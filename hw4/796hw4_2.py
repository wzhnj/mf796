# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 20:35:02 2021
Email:zuhuwang@bu.edu
@author: Zuhua Wang
"""

import pandas as pd

df = pd.read_excel('mf796-hw4-opt-data.xlsx')

print(df)

def a():
    #whether bid < ask
    for i in range(len(df)):
        if (df.iloc[i,3]>df.iloc[i,4]):
            print(df.iloc[i,])
        if (df.iloc[i,5]>df.iloc[i,6]):
            print(df.iloc[i,])
    #for the same K, longer time should have higher price
    for k in range(210,305,5):
        temp=df.loc[df['K']==k]
        temp_value=temp.iloc[0,3]

        for i in temp['call_bid']:
            if (i<temp_value):
                print("call_bid",k)
            temp_value=i
    for k in range(210,305,5):
        temp=df.loc[df['K']==k]
        temp_value=temp.iloc[0,4]

        for i in temp['call_ask']:
            if (i<temp_value):
                print("call_ask",k)
            temp_value=i
    for k in range(210,305,5):
        temp=df.loc[df['K']==k]
        temp_value=temp.iloc[0,5]

        for i in temp['put_bid']:
            if (i<temp_value):
                print("put_bid",k)
            temp_value=i
    for k in range(210,305,5):
        temp=df.loc[df['K']==k]
        temp_value=temp.iloc[0,6]

        for i in temp['put_ask']:
            if (i<temp_value):
                print("put_ask",k)
            temp_value=i
    #for the same expDays, K bigger, call price lower and put price higher
    for expDays in [49,140,203]:
        temp=df.loc[df['expDays']==expDays]
        temp_value=temp.iloc[0,3]

        for i in temp['call_bid']:
            if (i>temp_value):
                print("call_bid",k)
            temp_value=i
    for expDays in [49,140,203]:
        temp=df.loc[df['expDays']==expDays]
        temp_value=temp.iloc[0,4]

        for i in temp['call_ask']:
            if (i>temp_value):
                print("call_ask",k)
            temp_value=i
    for expDays in [49,140,203]:
        temp=df.loc[df['expDays']==expDays]
        temp_value=temp.iloc[0,5]

        for i in temp['put_bid']:
            if (i<temp_value):
                print("put_bid",k)
            temp_value=i
    for expDays in [49,140,203]:
        temp=df.loc[df['expDays']==expDays]
        temp_value=temp.iloc[0,6]

        for i in temp['put_ask']:
            if (i<temp_value):
                print("put_ask",k)
            temp_value=i
    #check convex
    for expDays in [49,140,203]:
        temp=df.loc[df['expDays']==expDays]
        for i in range(len(temp)-2):
            if (temp.iloc[i,3]+temp.iloc[i+2,3]-temp.iloc[i+1,3]<0):
                print("convex",expDays)
        for i in range(len(temp)-2):
            if (temp.iloc[i,4]+temp.iloc[i+2,4]-temp.iloc[i+1,4]<0):
                print("convex",expDays)
        for i in range(len(temp)-2):
            if (temp.iloc[i,5]+temp.iloc[i+2,5]-temp.iloc[i+1,5]<0):
                print("convex",expDays)
        for i in range(len(temp)-2):
            if (temp.iloc[i,6]+temp.iloc[i+2,6]-temp.iloc[i+1,6]<0):
                print("convex",expDays)
            
a()

































 