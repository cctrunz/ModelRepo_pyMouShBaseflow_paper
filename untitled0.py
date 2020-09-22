#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 09:32:28 2020

@author: cctrunz
"""

import numpy as np

stitch_initial = 122
stitch_total = 122 +100*6
addition = 0


for i in np.arange(70):
    if i % 2 == 0: #even The % sign is like division only it checks for the remainder, so if the number divided by 2 has a remainder of 0 it's even otherwise odd
        addition =  addition + 8          
    stitch_total = stitch_total + stitch_initial + addition
    if i % 5 == 0:
        stitch_total = stitch_total + stitch_initial + addition
        
    print(stitch_total)
