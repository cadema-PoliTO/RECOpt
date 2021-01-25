# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 23:40:09 2020

@author: giamm
"""

from time import time

##############################################################################

# The script is a Python adaption of Matlab's tic() toc().

_tstart_stack = []

def tic():
    _tstart_stack.append(time())
    
def toc():
    return(time() - _tstart_stack.pop())

# def tot_toc():
#     return(time() - _tstart_stack[0])


    
