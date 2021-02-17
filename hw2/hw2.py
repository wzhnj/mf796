# -*- coding: utf-8 -*-
"""
Created on Sun Jan 31 19:15:57 2021
Email:zuhuwang@bu.edu
@author: Zuhua Wang
"""
from math import *
from random import *
import numpy as np
from numpy import *

import scipy.stats as si
S0=10
r=0.01
sigma=0.2
T=3/12
K=12
def black_scholes_call(S0, K, T, r, sigma):
    ''' return put price calculated by black scholes model
    '''
    
    d1 = (np.log(S0 / K) + (r  + 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T))
    d2 = (np.log(S0 / K) + (r  - 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T))
    
    call = (S0 * si.norm.cdf(d1, 0.0, 1.0)-K * np.exp(-r * T) * si.norm.cdf(d2, 0.0, 1.0))
    #print("d1:",round(d1,4),"d2:",round(d2,4))
    return call

def left_riemann(nodes,x):
    if (x<-3):
        print("Error")
    dx=(x+3)/nodes
    result=si.norm.cdf(-3,0,1)
    start=(-3)
    for i in range(nodes):
        result=result+dx*si.norm.pdf(start)
        start=start+dx
    return result

def left_riemann_range(nodes,a,b):
    dx=(b-a)/nodes
    result=0
    start=a
    for i in range(nodes):
        result=result+dx*si.norm.pdf(start)
        start=start+dx
    return result

def midpoint_rule(nodes,x):
    if (x<-3):
        print("Error")
    dx=(x+3)/nodes
    result=si.norm.cdf(-3,0,1)
    start=(-3)+dx/2
    for i in range(nodes):
        result=result+dx*si.norm.pdf(start)
        start=start+dx
    return result

def midpoint_rule_range(nodes,a,b):
    
    dx=(b-a)/nodes
    result=0
    start=a+dx/2
    for i in range(nodes):
        result=result+dx*si.norm.pdf(start)
        start=start+dx
    return result

def calculation_error(x,result):
    return abs(si.norm.cdf(x)-result)

def bs_model_left_riemann(S0, K, T, r, sigma,nodes):
    d1 = (np.log(S0 / K) + (r  + 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T))
    d2 = (np.log(S0 / K) + (r  - 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T))
    call = (S0 * left_riemann(nodes,d1)-K * np.exp(-r * T) * left_riemann(nodes,d2))
    return call
def bs_model_midpoint_rule(S0, K, T, r, sigma,nodes):
    d1 = (np.log(S0 / K) + (r  + 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T))
    d2 = (np.log(S0 / K) + (r  - 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T))
    call = (S0 * midpoint_rule(nodes,d1)-K * np.exp(-r * T) * midpoint_rule(nodes,d2))
    #print(midpoint_rule(nodes,d1),midpoint_rule(nodes,d2))
    return call
def bs_model_gauss(S0, K, T, r, sigma,nodes):
    d1 = (np.log(S0 / K) + (r  + 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T))
    d2 = (np.log(S0 / K) + (r  - 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T))
    call = (S0 * (si.norm.cdf(-3)+gauss.GaussLegendreQuadrature(gauss.func , nodes, -3,d1)[0])-K * np.exp(-r * T) * (si.norm.cdf(-3)+gauss.GaussLegendreQuadrature(gauss.func , nodes, -3,d2)[0]))
    #print(si.norm.pdf(-3)+gauss.GaussLegendreQuadrature(gauss.func , nodes, -3,d1)[0],si.norm.pdf(-3)+gauss.GaussLegendreQuadrature(gauss.func , nodes, -3,d2)[0])
    return call

class gauss:
    ##################################################################
    # Recursive generation of the Legendre polynomial of order n
    def Legendre(n,x):
    	x=array(x)
    	if (n==0):
    		return x*0+1.0
    	elif (n==1):
    		return x
    	else:
    		return ((2.0*n-1.0)*x*gauss.Legendre(n-1,x)-(n-1)*gauss.Legendre(n-2,x))/n
        
        
    def DLegendre(n,x):
    	x=array(x)
    	if (n==0):
    		return x*0
    	elif (n==1):
    		return x*0+1.0
    	else:
    		return (n/(x**2-1.0))*(x*gauss.Legendre(n,x)-gauss.Legendre(n-1,x))
    
    
    def LegendreRoots(polyorder, tolerance=1e-20):
        if (polyorder<2):
            err=1
        else:
            roots=[]
            for i in range(1, int(polyorder/2)+1):
                x=cos(pi*(i-0.25)/(polyorder+0.5))
                error=10*tolerance
                iters=0
                while(error>tolerance) and (iters<1000):
                    dx=-gauss.Legendre(polyorder,x)/gauss.DLegendre(polyorder,x)
                    x=x+dx
                    iters=iters+1
                    error=abs(dx)
                roots.append(x)
            roots= array(roots)
            if (polyorder%2==0):
                roots=concatenate((-1.0*roots,roots[::-1]))
            else:
                roots=concatenate((-1.0*roots, [0.0], roots[::-1]))
            err=0
        return[roots,err]
        
        
    def GaussLegendreWeights(polyorder):
    	W=[]
    	[xis,err]=gauss.LegendreRoots(polyorder)
    	if err==0:
    		W=2.0/( (1.0-xis**2)*(gauss.DLegendre(polyorder,xis)**2) )
    		err=0
    	else:
    		err=1 # could not determine roots - so no weights
    	return [W, xis, err]
    # The integral value 
    # func 		: the integrand
    # a, b 		: lower and upper limits of the integral
    # polyorder 	: order of the Legendre polynomial to be used
    #
    def GaussLegendreQuadrature(func, polyorder, a, b):
    	[Ws,xs, err]= gauss.GaussLegendreWeights(polyorder)
    	if (err==0):
    		ans=(b-a)*0.5*sum( Ws*func( (b-a)*0.5*xs+ (b+a)*0.5 ) )
    	else: 
    		# (in case of error)
    		err=1
    		ans=None
    	return [ans,err]
        
    # The integrand - change as required
    def func(x):
        #return 1
    	return exp(-0.5*((x-0.0)/1)**2)/(1*sqrt(2*pi))

def estimate_error():
    start=0
    start1=-1
    start2=-2
    end1=1
    end2=2
    left1=abs((si.norm.cdf(end1)-si.norm.cdf(start))-(left_riemann_range(500,start,end1)))
    left2=abs((si.norm.cdf(end2)-si.norm.cdf(start))-(left_riemann_range(500,start,end2)))
    print("left riemann length 1 error:",left1)
    print("left riemann length 2 error:",left2)
    left50=abs((si.norm.cdf(end1)-si.norm.cdf(start))-(left_riemann_range(50,start,end1)))
    left500=abs((si.norm.cdf(end1)-si.norm.cdf(start))-(left_riemann_range(500,start,end1)))
    print("left riemann 50 dx error:",left50)
    print("left riemann 500 dx error",left500)
    
    mid1=abs((si.norm.cdf(end1)-si.norm.cdf(start1))-(midpoint_rule_range(500,start1,end1)))
    mid2=abs((si.norm.cdf(end2)-si.norm.cdf(start2))-(midpoint_rule_range(500,start2,end2)))
    print("midpoint rule length 1 error:",mid1)
    print("midpoint rule length 2 error:",mid2)
    mid50=abs((si.norm.cdf(end1)-si.norm.cdf(start))-(midpoint_rule_range(50,start,end1)))
    mid500=abs((si.norm.cdf(end1)-si.norm.cdf(start))-(midpoint_rule_range(500,start,end1)))
    print("midpoint rule 50 dx error:",mid50)
    print("midpoint rule 500 dx error",mid500)
    print(mid2/mid1)
    
    gauss1=abs((si.norm.cdf(end1)-si.norm.cdf(start1))-(gauss.GaussLegendreQuadrature(gauss.func , 5, start,end1)[0]))
    gauss2=abs((si.norm.cdf(end2)-si.norm.cdf(start2))-(gauss.GaussLegendreQuadrature(gauss.func , 5, start,end2)[0]))
    print("GaussLegendreQuuadrature length 1 error:",gauss1)
    print("GaussLegendreQuuadrature length 2 error:",gauss2)
    gauss3=abs((si.norm.cdf(end1)-si.norm.cdf(start))-(gauss.GaussLegendreQuadrature(gauss.func , 3, start,end1)[0]))
    gauss5=abs((si.norm.cdf(end1)-si.norm.cdf(start))-(gauss.GaussLegendreQuadrature(gauss.func , 5, start,end1)[0]))
    print("GaussLegendreQuuadrature 3 dx error:",gauss3)
    print("GaussLegendreQuuadrature 5 dx error",gauss5)
    print(gauss3/gauss5)
    
#estimate_error()
'''
print("Problem 1.1: The price is",black_scholes_call(S0,K,T,r,sigma))
print("Problem 1.2.a: For N=5, the error is",calculation_error(0,left_riemann(5,0)))
print("For N=10, the error is",calculation_error(0,left_riemann(10,0)))
print("For N=50, the error is",calculation_error(0,left_riemann(50,0)))
print("For N=5, the price is",bs_model_left_riemann(S0, K, T, r, sigma,5))
print("For N=10, the price is",bs_model_left_riemann(S0, K, T, r, sigma,10))
print("For N=50, the price is",bs_model_left_riemann(S0, K, T, r, sigma,50))
print("Problem 1.2.b: For N=5, the error is",calculation_error(0,midpoint_rule(5,0)))
print("For N=10, the error is",calculation_error(0,midpoint_rule(10,0)))
print("For N=50, the error is",calculation_error(0,midpoint_rule(50,0)))
print("For N=5, the price is",bs_model_midpoint_rule(S0, K, T, r, sigma,5))
print("For N=10, the price is",bs_model_midpoint_rule(S0, K, T, r, sigma,10))
print("For N=50, the price is",bs_model_midpoint_rule(S0, K, T, r, sigma,50))
print("Problem 1.2.c: For N=5, the error is",calculation_error(0,si.norm.cdf(-3)+gauss.GaussLegendreQuadrature(gauss.func , 5, -3,0)[0]))
print("For N=10, the error is",calculation_error(0,si.norm.cdf(-3)+gauss.GaussLegendreQuadrature(gauss.func , 10, -3,0)[0]))
print("For N=15, the error is",calculation_error(0,si.norm.cdf(-3)+gauss.GaussLegendreQuadrature(gauss.func , 15, -3,0)[0]))
print("For N=5, the price is",bs_model_gauss(S0, K, T, r, sigma,5))
print("For N=10, the price is",bs_model_gauss(S0, K, T, r, sigma,10))
print("For N=15, the price is",bs_model_gauss(S0, K, T, r, sigma,15))
'''
'''
def test():
    ratio=[]
    for i in range(100):
        start=uniform(-10,10)
        length=uniform(0,10)
        end1=start+length
        end2=start+length*2
        left1=abs((si.norm.cdf(end2)-si.norm.cdf(start))-(left_riemann_range(500,start,end2)))
        left2=abs((si.norm.cdf(end1)-si.norm.cdf(start))-(left_riemann_range(500,start,end1)))
        ratio.append(left2/left1)
    print(median(ratio))
test()
'''
'''
left1=abs((si.norm.cdf(-2)-si.norm.cdf(-3))-(midpoint_rule_range(500,-3,-2)))
left2=abs((si.norm.cdf(-1)-si.norm.cdf(-3))-(midpoint_rule_range(500,-3,-1)))
print("ratio:",left2/left1)
left1=abs((si.norm.cdf(-2)-si.norm.cdf(-3))-gauss.GaussLegendreQuadrature(gauss.func,3,-3,-2)[0])
left2=abs((si.norm.cdf(-2)-si.norm.cdf(-3))-gauss.GaussLegendreQuadrature(gauss.func,6,-3,-2)[0])
print(left2,left1)
'''




'''
print(calculation_error(-1,left_riemann(5,-1)))
print(calculation_error(-1,left_riemann(5,-1))/calculation_error(-2,left_riemann(5,-2)))
print(calculation_error(-1,midpoint_rule(10,-1)))
print(calculation_error(-2,midpoint_rule(10,-2)))
print(calculation_error(-1,midpoint_rule(10,-1))/calculation_error(-2,midpoint_rule(10,-2)))

'''
#print(3**3)
#print(calculation_error(0,midpoint_rule(10,0)))
#print(calculation_error(0,si.norm.cdf(-3)+gauss.GaussLegendreQuadrature(gauss.func , 10, -3,0)[0]))

def midpoint_rule_range_inner(nodes,a,b,s2):
    
    dx=(b-a)/nodes
    result=0
    start=a+dx/2
    for i in range(nodes):
        result=result+dx*(start-380)/(2*pi*20*15*sqrt(1-0.95**2))*exp(-1/(2*(1-0.95**2))*(((start-380)/20)**2+((s2-380)/15)**2-(2*0.95*(start-380)*(s2-380)/(20*15))))
        start=start+dx
    return result

def midpoint_rule_range_outer(nodes,a,b,c,d):
    dx=(b-a)/nodes
    result=0
    start=a+dx/2
    for i in range(nodes):
        result=result+dx*midpoint_rule_range_inner(nodes,c,d,start)
        start=start+dx
    return result

def b_b():
    #print(midpoint_rule_range_outer(10,375,500,380,500))
    #print(midpoint_rule_range_outer(20,375,500,380,500))
    #print(midpoint_rule_range_outer(400,375,900,380,900))
    print(midpoint_rule_range_outer(4000,0,360,380,1000))



b_b()






