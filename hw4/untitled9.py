# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 13:52:15 2021
Email:zuhuwang@bu.edu
@author: Zuhua Wang
"""


from pandas_datareader import data 
from datetime import date
import yfinance as yf
yf.pdr_override()
import pandas as pd
tickers = ['AAPL', 'MSFT', '^GSPC']
import csv





# We would like all available data from 01/01/2000 until 12/31/2016.
start_date = '2011-01-01'
end_date = '2015-12-31'
df=pd.DataFrame()
# User pandas_reader.data.DataReader to load the desired data. As simple as that.
for i in range(len(tickers)):
    df[tickers[i]] = data.DataReader(tickers[i], start_date, end_date).iloc[:,3]
    
    
print(df)

