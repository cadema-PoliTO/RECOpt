# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 12:12:00 2020

@author: giamm
"""

import numpy as np

###############################################################################


# This script contains a method that receives as input one or more load profiles
# (a column vector or morcolumn vector horizotally stacked in a 2d-array) for a 
# 1440 min simulation time, with a resolution (timestep) of 1 min. It returns
# the load profile with a different resolution, that is given as input (from
# 5 min to 60 min), keeping constant the total energy consumed in each time
# interval.


###############################################################################

def aggregator(load_profile,dt_aggr):
    
    '''The function aggregates a load_profile, changing the time resolution.
    
    Inputs:
        load_profile - load-profile to be aggregated, with a resolution of 1 min (W)
        dt_aggr - new timestep (time resolution) to be used (min)
        
    Output:
        lp_aggr - aggregated load profile (W)
    '''
    
    # In case load_profile is a 1d array, it is reshaped into a 2d array (as a
    # column vector) so that the function can work properly.
    
    flag = len(np.shape(load_profile))
    if flag == 1:
        load_profile = load_profile.reshape((np.size(load_profile,axis=0),1))
   
    ##########  Original time discretization
    
    dt = 1 #original timestep (min)
    time = 1440 #total simulation time (min)
    
    
    ########## New time discretization
    
    length = np.size(load_profile,axis=0) #length of each columns vector
    length_aggr = int(length*dt/dt_aggr) #length of the new load profile vector
    
    # The new time-step needs to be a submultiple of the total time (1440 min)
    # other-wise the last time interval will not end exactly at the end of 23:59
    # If thid does not happen, the function returns an error message.
    if length_aggr != length*dt/dt_aggr:
        print('Choose a different aggregation timestep, please')
        return
    
    # The load-profile is aggregated with a different timestep. Since both the
    # original vector and the new one are evenly distributed on time (meaning 
    # that the timestep is the same for all the time intervals), this is easily 
    # done with only one "for" loop and with no need to know the time vectors.
    # First, an array is created where to store the load profile(s) with the 
    # new time discretization.     
    lp_aggr = np.zeros((length_aggr,np.size(load_profile,axis=1)))
    jj = 0
    
    for ii in range(0,length_aggr):
        
        lp_aggr[ii,:] = sum(load_profile[jj:jj + int(dt_aggr/dt),:])*(dt/dt_aggr)
        jj = jj + int(dt_aggr/dt)
      
        
    # In case load_profile was a 1d array, the new one is reshaped into a 
    # 2d array (as a column vector) so that a vector is returned
    if flag == 1:
        lp_aggr = lp_aggr.reshape((np.size(lp_aggr,axis=0),))
        
    return(lp_aggr)


# # Uncomment the following lines to test the function
# import matplotlib.pyplot as plt 
# dt = 1
# time = np.arange(0,100,dt)
# powers = np.random.randint(100,size=(100,1))

# dt_aggr = 10
# time_aggr = np.arange(0,100,10)
# powers_aggr = aggregator(powers,dt_aggr)

# for ii in range(0,np.size(powers,axis=1)):
    
#     plt.figure()
#     plt.bar(time,powers[:,ii],width=dt,align='edge')
#     plt.bar(time_aggr,powers_aggr[:,ii],width=dt_aggr,align='edge',fill=False,edgecolor='k')
    
#     print(np.sum(powers[:,ii])*dt)
#     print(np.sum(powers_aggr[:,ii])*dt_aggr)
    
    