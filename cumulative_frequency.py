# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 12:09:59 2020

@author: giamm
"""

import numpy as np

##############################################################################


# This script is used to create a function that receives as inputs a time-vector
# (representing the x-axis) and a frequency density vector (representing the y-axis)
# and returns the cumulative frequency.

##############################################################################

def cum_freq(time,freqdens):
      
    '''The function computes the cumulative frequency.
    
    Inputs:
        time - time vector on x-axis (min)
        freqdens - frequency density on y-axis (W/min)

        
    Output:
        cumfreq - cumulative frequency
    '''
    
    # Creating the cumulate by adding the single values of the frequency density in each time interval 
    res = np.zeros(np.shape(time))
    cum = 0
    for ii in range(1,np.size(time)):
                
        res[ii] = 0.5*(freqdens[ii] + freqdens[ii-1])*(time[ii] - time[ii-1]) + cum
        cum = res[ii]
            
    # Normalization and return
    return(res/cum)  

    
# ##### Function testing: uncomment the following line to test the function
 
# import matplotlib.pyplot as plt
# from profile_interpolation import interp_profile
  
# dt = 10
# dt_intp = 1
# time = np.arange(0,100,dt) #(h)
# power = np.random.randint(100,size=np.size(time)) #(W)
# time_intp = np.arange(0,100,dt_intp)
# power_intp = interp_profile(time,power,time_intp)

# plt.figure()
# plt.bar(time,power,width=dt,align='edge')
# plt.plot(time_intp,power_intp,'r')

# res = cum_freq(time,power)
# res_intp = cum_freq(time_intp,power_intp)

# res1 = np.cumsum(power)/np.sum(power)
# res_intp1 = np.cumsum(power_intp)/np.sum(power_intp)

# res_intp2 = np.zeros(np.shape(time_intp));
# cum = 0;
# for ii in range(1,np.size(time_intp)):
                
#     res_intp2[ii] = 0.5*(power_intp[ii] + power_intp[ii-1])*(time_intp[ii] - time_intp[ii-1]) + cum
#     cum = res_intp2[ii]  
# res_intp2 = res_intp2/cum

# plt.figure()
# plt.bar(time,res,width=dt,align='edge')
# plt.bar(time,res1,width=dt,align='edge', fill = False, edgecolor='r')
# plt.plot(time_intp,res_intp,'b')
# plt.plot(time_intp,res_intp,'r')
# plt.plot(time_intp,res_intp2,'k--')



    
    