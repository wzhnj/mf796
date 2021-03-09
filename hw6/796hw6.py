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
from numpy import linalg as LA


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

#max price
x_max=500
x_min=0
#difference price
Nspace=200
#difference time
Ntime=3000
#risk free rate
r=0.05/100
#strike price
K=390
start = dt.date( 2021, 3, 3 )
end = dt.date( 2021, 9, 30 )
days = np.busday_count( start, end )
#total time 
Texpir=days/252

x, dx = np.linspace(x_min, x_max, Nspace, retstep=True)   # space discretization
T, dt = np.linspace(0, Texpir, Ntime, retstep=True)       # time discretization
Payoff = np.maximum(x-K,0)          # Call payoff

V = np.zeros((Nspace,Ntime))       # grid initialization
offset = np.zeros(Nspace-2)        # vector to be used for the boundary terms   

V[:,-1] = Payoff                   # terminal conditions 

sig2 = sigma**2; dxx = dx**2


a=1-sig2*(x**2)*dt/dxx-r*dt
l=sig2*(x**2)*dt/2/dxx-r*x*dt/2/dx
u=sig2*(x**2)*dt/2/dxx+r*x*dt/2/dx

A=np.zeros((Nspace-2,Nspace-2))

A[0][0]=a[1]
A[0][1]=u[1]
A[-1][-2]=l[-2]
A[-1][-1]=a[-2]
#M=Matrix(A)
print(A)


for i in range(1,Nspace-3):
    A[i][i-1]=l[i+1]
    #print(i,'a')
    A[i][i]=a[i+1]
    #print(i,'b')
    A[i][i+1]=u[i+1]
    #print(i,'c')
print("A");print(A)

print('V');print(V)
V=np.matrix(V)
count=1
c=0
for i in range(Ntime-2,-1,-1):
    b=u[-2]*(x[-1]-K*math.exp(-1*r*dt*count))
    count+=1
    V[1:-1,i]=np.matrix(A)*np.matrix(V[1:-1,i+1])
    V[-2,i]+=b
    #if (i<4):
     #   break
    
    for j in range(1,Nspace):
        
        if (V[j,i]-(x[j]-K)<0):
            #print(1)
            c=c+1
        else:
            #print(2)
        
        V[j,i]=max(V[j,i],x[j]-K)
'''
           print(j)
           print(V[j,i]-(x[j]-K))
           #print()
           print(V[j,i],x[j]-K)
'''     
print(c)
'''
print(V[:,[0,1,2]])
i=-4
print(V[:,[i-2,i-1,i]])
    
print(abs(sig2*x_max**2*dt/dx**2))
#median=int(Nspace/2)
print(V[25,0])
#print(median)
    
eigenvalue, _ = LA.eig(A)
print("max(eigen_value)",max(eigenvalue))   
'''
#for i in range(len(x)):
#    print(i, x[i],V[i,0])
#print(eigenvalue)

