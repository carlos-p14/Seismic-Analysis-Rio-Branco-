# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 19:33:27 2021

@author: Carlos

This script will only interpolate the Fragility Functions
"""

import numpy as np

def interpolate(x, y, x0):
    
    if x0 > max(x):
        return(max(y))
        
    elif x0 < min(x):
        return(np.interp(x0, [0, min(x)], [0, min(y)]))
    
    else:
        return(np.interp(x0,x,y))    