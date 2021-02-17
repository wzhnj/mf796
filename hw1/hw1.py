# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 22:34:48 2021
Email:zuhuwang@bu.edu
@author: Zuhua Wang
"""
import numpy as np
import scipy.stats as si
S0=100
r=0.0
beta=1
sigma=0.4
T=1
K=100
def way():
    way=[S0]
    for i in range(252):
        way.append(way[i]+way[i]*r*(1/252)+(way[i]**beta)*sigma*np.random.normal(0,(1/252)**0.5))
    return way[252]
def b():
    record=[]
    temp=0
    for i in range(100):
        temp=way()
        if (temp>=100):
            record.append(temp-100)
        else:
            record.append(0)
    print(record)
    print(np.mean(record))
    return np.mean(record)
def black_scholes_call(S0, K, T, r, sigma):
    ''' return put price calculated by black scholes model
    '''
    
    d1 = (np.log(S0 / K) + (r  + 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T))
    d2 = (np.log(S0 / K) + (r  - 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T))
    
    call = (S0 * si.norm.cdf(d1, 0.0, 1.0)-K * np.exp(-r * T) * si.norm.cdf(d2, 0.0, 1.0))
    #print("delta= ",si.norm.cdf(d1, 0.0, 1.0))
    return call
def delta(S0, K, T, r, sigma):
    d1 = (np.log(S0 / K) + (r  + 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T))
    return si.norm.cdf(d1, 0.0, 1.0)
def c():
    #print(black_scholes_call(S0, K, T, r, sigma))
    return black_scholes_call(S0, K, T, r, sigma)
def f():
    payofflist=[]
    for i in range(10000):
        price=way()
        if (price>=100):
            payoff=price-100
        else:
            payoff=0
        payoff=payoff-(price-100)
        payofflist.append(payoff)
    #print(payofflist)
    #print(np.mean(payofflist))
    return np.mean(payofflist)
c()
print(f())
'''
l=[]
for i in range(100):
    l.append(f())
print(np.mean(l))
'''