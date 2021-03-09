# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 13:52:15 2021
Email:zuhuwang@bu.edu
@author: Zuhua Wang
"""
from sklearn.decomposition import PCA
from datapackage import Package
from pandas_datareader import data 
from datetime import date
import yfinance as yf
yf.pdr_override()
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
from numpy.linalg import inv

'''
tickers = ['AAPL', 'MSFT', '^GSPC']
import csv
with open('tickers.csv') as f:
    reader=csv.reader(f)
    tickers = list(reader)
tickers=tickers[0:1]
#print(tickers)
'''


package = Package('https://datahub.io/core/s-and-p-500-companies/datapackage.json')

# print list of all resources:
#print(package.resource_names)
tickers=[]
# print processed tabular data (if exists any)
for resource in package.resources:
    if resource.descriptor['datahub']['type'] == 'derived/csv':
        tickers=resource.read()


# We would like all available data from 01/01/2000 until 12/31/2016.
start_date = '2015-01-01'
end_date = '2019-12-31'
df=pd.DataFrame()
# User pandas_reader.data.DataReader to load the desired data. As simple as that.
for i in range(104):
    if (tickers[i][0]=='CARR' or tickers[i][0]=='BF.B' or tickers[i][0]=='BRK.B' or tickers[i][0]=='ABC'):
        continue
    df[tickers[i][0]] = data.DataReader(tickers[i][0], start_date, end_date).iloc[:,3]
    
    
print(df)
df.isnull().values.any()
#2
df=df.pct_change(1)
df=df.iloc[1:len(df),:]
cov=df.cov()
print(np.linalg.eigvals(cov))
#4
pca=PCA()
pca.fit(df)
print(pca.explained_variance_ratio_)

explained_ratio=pca.explained_variance_ratio_
ratio=0
for i in range(len(explained_ratio)):
    ratio=ratio+explained_ratio[i]
    print(i, ratio)
pca=PCA(55)
pca.fit(df)
scores=pca.transform(df)
reconstruct = pca.inverse_transform(scores)
residual=df-reconstruct
print(residual)



plt.figure()
residual.plot()





r1=[]
r2=[]
for i in range(100):
    r1.append(1)
    if (i<17):
        r2.append(1)
    else:
        r2.append(0)

g=np.matrix([r1, r2])

c_inv=inv(cov)

g_tran = np.transpose(g)

gcg = np.dot(np.dot(g,c_inv),g_tran)

c=np.matrix([[1],[0.1]])

r=df.mean()
lam=np.dot(inv(gcg),2*c-np.transpose(np.dot(np.dot(g,c_inv),r)))

a=1

w = np.dot(1/(2*a)*c_inv,(np.transpose(np.matrix(r))-np.dot(g_tran,lam)))

print(w)










