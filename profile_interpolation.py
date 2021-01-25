# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 16:16:02 2020

@author: giamm
"""

import numpy as np

##############################################################################


# This script is used to create a profile, interpolating a set of data


##############################################################################
    
def interp_profile(time,profile,time_intp):
      
    '''The returns a linear-interpolated profile (with periodic behavior).
    
    Inputs:
        time - time vector (h)
        profile - profile to be interpolated
        time_intp - new time vector (h)
        
    Output:
        profile_intp - interpolated profile
    '''
       
    ########## Linear interpolation

    dt = (time[-1]-time[0])/(np.size(time)-1) #timestep of time vector
    period = time[-1] + dt #period 
    
    #Linear interpolation
    profile_intp  = np.interp(time_intp, time, profile, period = period)
        
    # Saturation to 0 for negative values (since negative powers don't have physical meaning)
    profile_intp[profile_intp < 0] = 0
    
    return(profile_intp)


# # Function testing: uncomment the following line to test the function
# from scipy.interpolate import CubicSpline
# from scipy.interpolate import interp1d
# import datareader 
# import matplotlib.pyplot as plt
# filename = 'freqdens_wm_we_sawp.csv'      
# data = datareader.read_general(filename,';')
# time = data[:,0] #(h)
# time = time*60 #conversion to minutes
# power = data[:,1] #(W)
# time_intp = np.arange(0,1440,1)

# # time = np.arange(0,100,10)
# # power = np.random.randint(100,size = np.size(time))
# # time_intp = np.arange(0,100,1)

# dt = (time[-1]-time[0])/(np.size(time)-1)
# dt_intp = (time_intp[-1]-time_intp[0])/(np.size(time_intp)-1)

# power_intp = interp_profile(time,power,time_intp)

# f = CubicSpline(np.append(time,time[-1]+(time[-1]-time[-2])),np.append(power,power[0]),bc_type='periodic')
# power_intp1 = f(time_intp)

# f = interp1d(time,power,fill_value='extrapolate')
# power_intp2 = f(time_intp)

# power_intp3 = np.interp(time_intp,time,power,period = 1440)

# plt.figure()
# plt.bar(time,power,align='edge',width=dt)
# plt.plot(time_intp,power_intp1,'b')
# plt.plot(time_intp,power_intp2,'r')
# plt.plot(time_intp,power_intp3,'k')

# plt.figure()
# plt.plot(time_intp,power_intp1,'b')
# plt.plot(time_intp,power_intp2,'r')
# plt.plot(time_intp,power_intp3,'k')

# plt.xlim([0,100])

# plt.figure()
# plt.plot(time_intp,power_intp1,'b')
# plt.plot(time_intp,power_intp2,'r')
# plt.plot(time_intp,power_intp3,'k--')

# plt.xlim([1340,1439])


# en=sum(power)*(time[1]-time[0])/60 #Wh
# en_intp=np.trapz(power_intp1,time_intp)/60
# en_intp2=np.trapz(power_intp2,time_intp)/60
# en_intp3=np.trapz(power_intp3,time_intp)/60
