# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 17:24:47 2021
Email:zuhuwang@bu.edu
@author: Zuhua Wang
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 21:21:49 2021
Email:zuhuwang@bu.edu
@author: Zuhua Wang
"""

import yfinance as yf
import math
import numpy as np
import datetime as dt

SPY = yf.Ticker("SPY")

# get stock info


# get historical market data
hist = SPY.history(period="1y")
data=hist['Close']
#print(data)
logret=[]
for i in data:
    logret.append(math.log(i))
#print(logret)
z=[]
for i in range(len(logret)-1):
    z.append(logret[i+1]-logret[i])
#print(z)
sigma = np.std(z)*252**0.5
#print()
print("sigma",sigma)

#problem 4

M=387
hs=1
ht=1/252
r=0.05/100
K=385
start = dt.date( 2021, 3, 3 )
end = dt.date( 2021, 9, 30 )

days = np.busday_count( start, end )

tend=days/252

price=np.zeros((M+1,days+1))

s=np.array(range(M+1))
ct=s-K

for i in range(len(ct)):
    if (ct[i]<0):
        ct[i]=0
price[:,days]=ct
price=np.matrix(price)

a=1-(sigma**2)*(s**2)*ht/(hs**2)-r*ht
l=(sigma**2)*(s**2)*ht/2/(hs**2)-r*s*ht/2/hs
u=(sigma**2)*(s**2)*ht/2/(hs**2)+r*s*ht/2/hs

A=np.zeros((M-1,M-1))

A[0][0]=a[1]
A[0][1]=u[1]
A[M-1-1][M-2-1]=l[M-1]
A[M-1-1][M-1-1]=a[M-1]



for i in range(1,M-2):
    A[i][i-1]=l[i+1]
    #print(i,'a')
    A[i][i]=a[i+1]
    #print(i,'b')
    A[i][i+1]=u[i+1]
    #print(i,'c')


print(price[:,-1])
for i in range(days):
    b=u[M-1]*(s[M]-K*math.exp(-1*r*((i+1)*(1/252))))
    price[1:M,days-(i+1)]=np.matrix(A)*np.matrix(price[1:M,days-i])
    price[M-1,days-(i+1)]+=b
    #if (i>3):
        #break
print(price[:,[-3,-2]])
    
w, v = LA.eig(A)
print(w)
print(v)
    

    
    
    
    
    
    
    