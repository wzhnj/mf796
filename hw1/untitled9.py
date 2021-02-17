# -*- coding: utf-8 -*-
"""
Created on Sun Jan 31 19:15:57 2021
Email:zuhuwang@bu.edu
@author: Zuhua Wang
"""

import smtplib, ssl
import time

port = 465  # For SSL
password = "@Quanlan1"
# Create a secure SSL context
context = ssl.create_default_context()
a=106
'''
for i in range(3):    
    with smtplib.SMTP_SSL("smtp.gmail.com", port, context=context) as server:
        server.login("zuhuawang2021@gmail.com", password)
        # TODO: Send email here
        sender_email = "zuhuawang2021@gmail.com"
        receiver_email = "zuhuwang@bu.edu"
        message = "There is only "+str(a)+" days to go home."
        server.sendmail(sender_email, receiver_email, message)
        print(a)
        a-=1
        time.sleep(10)
        '''
print(context)