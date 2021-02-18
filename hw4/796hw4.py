# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 20:46:30 2021
Email:zuhuwang@bu.edu
@author: Zuhua Wang
"""
import pandas as pd
import numpy as np
from scipy.stats import norm
from scipy import optimize
from scipy.optimize import fsolve
import math
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

'''A = np.array([[1,2],[3,4]])
B=np.array([1,2])
X = np.linalg.inv(A).dot(B)
print(X)
'''
df1=pd.DataFrame(np.array([[32.25,28.36],[24.73,21.78],
                           [20.21,18.18],[18.24,16.45],
                           [15.74,14.62],[13.70,12.56],
                           [11.48,10.94]]),columns=[1,2])
delta = [10/100,25/100,40/100,50/100,40/100,25/100,10/100]
put_call = ['p','p','p','c','c','c','c']
df1["put_call"]=put_call
df1['delta']=delta
def a():
    '''
    add strike column
    '''
    
    one_strike=[]
    three_strike=[]
    for i in range(len(df1)):
        '''
        find strike
        '''
        if (df1.iloc[i,2]=="p"):
            delta = -df1.iloc[i,3]
            nd1=delta + 1
            d1=norm.ppf(nd1)
            sigma=df1.iloc[i,0]/100
            def f(x):
                #return (3*x-1)
                return (d1-(math.log(100/x)+(1/12)*(0+(sigma**2)/2))/(sigma*math.sqrt(1/12)))
            root=fsolve(f, 100)
            #k=100/math.exp(d1*sigma*math.sqrt(1/12)-(1/12)*sigma*sigma/2)
            #print(nd1,d1,sigma,k)
            one_strike.append(root)
        else:
            delta = df1.iloc[i,3]
            nd1=delta
            d1=norm.ppf(nd1)
            sigma=df1.iloc[i,0]/100
            def f(x):
                #return (3*x-1)
                return (d1-(math.log(100/x)+(1/12)*(0+(sigma**2)/2))/(sigma*math.sqrt(1/12)))
            root=fsolve(f, 100)
            #k=100/math.exp(d1*sigma*math.sqrt(1/12)-(1/12)*sigma*sigma/2)
            #print(nd1,d1,sigma,k)
            one_strike.append(root)
    for i in range(len(df1)):
        '''
        find strike
        '''
        if (df1.iloc[i,2]=="p"):
            delta = -df1.iloc[i,3]
            nd1=delta + 1
            d1=norm.ppf(nd1)
            sigma=df1.iloc[i,1]/100
            def f(x):
                #return (3*x-1)
                return (d1-(math.log(100/x)+(1/12)*(0+(sigma**2)/2))/(sigma*math.sqrt(1/12)))
            root=fsolve(f, 100)
            #k=100/math.exp(d1*sigma*math.sqrt(1/12)-(1/12)*sigma*sigma/2)
            #print(nd1,d1,sigma,k)
            three_strike.append(root)
        else:
            delta = df1.iloc[i,3]
            nd1=delta
            d1=norm.ppf(nd1)
            sigma=df1.iloc[i,1]/100
            def f(x):
                #return (3*x-1)
                return (d1-(math.log(100/x)+(1/12)*(0+(sigma**2)/2))/(sigma*math.sqrt(1/12)))
            root=fsolve(f, 100)
            #k=100/math.exp(d1*sigma*math.sqrt(1/12)-(1/12)*sigma*sigma/2)
            #print(nd1,d1,sigma,k)
            three_strike.append(root)
        
    df1['three_strike']=three_strike
    df1['one_strike']=one_strike
    print (df1)
    return
def b():
    x = df1.iloc[:,5].to_numpy()
    y = df1.iloc[:,0].to_numpy()
    f = interp1d(x, y)
    print(type(x),y)
    xnew = np.linspace(89.2, 104, num=41, endpoint=True)
    plt.plot(x, y, 'o', xnew, f(xnew), '-')
    plt.legend(['data', 'linear'], loc='best')
    plt.title('1M')
    plt.xlabel('strike')
    plt.ylabel('volatility')
    plt.show()
def b3():
    x = df1.iloc[:,4].to_numpy()
    y = df1.iloc[:,1].to_numpy()
    f = interp1d(x, y)
    print(type(x),y)
    xnew = np.linspace(90.4, 104.17, num=41, endpoint=True)
    plt.plot(x, y, 'o', xnew, f(xnew), '-')
    plt.legend(['data', 'linear'], loc='best')
    plt.title('3M')
    plt.xlabel('strike')
    plt.ylabel('volatility')
    plt.show()
    
def call(S, K, T, r, sigma):
    
    #S: spot price
    #K: strike price
    #T: time to maturity
    #r: interest rate
    #sigma: volatility of underlying asset
    
    d1 = (np.log(S / K) + (r + 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T))
    d2 = (np.log(S / K) + (r - 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T))
    
    call = (S * si.norm.cdf(d1, 0.0, 1.0) - K * np.exp(-r * T) * si.norm.cdf(d2, 0.0, 1.0))
    
    return call
def c():
    k=list(range(1, 150))
    phi = (call()-call()*2+call())/(0.1)**2

a()
b()
b3()
#print(df1)