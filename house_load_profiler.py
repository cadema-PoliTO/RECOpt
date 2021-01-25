# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 18:37:58 2020

@author: giamm
"""

# from tictoc import tic, toc

import numpy as np

import datareader #routine created to properly read the files needed in the following 
from load_profiler import load_profiler as lp

###############################################################################


# This file contains a method that, for a given household (considering its 
# availability of each appliance considered in this study) returns the electric
# load profile for the household during one day (1440 min), with a resolution
# of 1 minute.


###############################################################################

########## Routine

def house_load_profiler(apps_availability, day, season, appliances_data, **params):
    
    ''' The method returns a load profile for a given household in a total simulation time of 1440 min, with a timestep of 1 min.
    
    Inputs:
        apps_availability - 1d-array, availability of each appliance for the household is stored (1 if present, 0 if not)
        day - str, type of day (weekday: 'wd'| weekend: 'we')
        season - str, season (summer: 's', winter: 'w', autumn or spring: 'ap')
        appliances_data - dict, various input data related to the appliances
        params - dict, simulation parameters
        
    Outputs:
        house_load_profile - 1d-array, load profile for the household (W)
        energy - 1d-array, energy consumed by each appliance in one day (Wh/day)
        
    '''

    ## Time 
    # Time discretization for the simulation    

    # Time-step, total time and vector of time from 00:00 to 23:59 (one day) (min)
    dt = 1 
    time = 1440 
    time_sim = np.arange(0,time,dt) 

    
    ## Parameters
    # Simulation parameters that can be changed from the user

    # Energy class of the appliances
    power_max = params['power_max']

   
    ## Input data for the appliances
    # Appliances' attributes, energy consumptions and user's coefficients 

    # apps is a 2d-array in which, for each appliance (rows) and attribute value is given (columns)
    apps_ID = appliances_data['apps_ID']

    # apps_attr is a dictionary in which the name of each attribute (value) is linked to its columns number in apps (key)
    apps_attr = appliances_data['apps_attr']


    ## Household's load profile
    # Generating the load profile for the house, considering which appliances are available

    # Initializing the power vector for the load profile (W)
    house_load_profile = np.zeros(np.shape(time_sim)) 

    # Initializing the vector where to store the energy consumption from each appliance (Wh/day)
    energy = np.zeros(len(apps_ID)) 
    
    # Using the method load_profiler(lp) to get the load profile for each appliance    
    for app in apps_ID:  

        # The ID number of the appliance is stored in a variable since it will be used man times
        app_ID = apps_ID[app][apps_attr['id_number']]     
        
        # Skipping appliances that are not present in the household
        if apps_availability[app_ID] == 0:
            continue
        
        load_profile = lp(app, day, season, appliances_data, **params) #load_profile has to outputs (time and power)
        
        # In case the instantaneous power exceedes the maximum power, some tries
        # are made in order to change the moment in which the next appliance is
        # switched on (since it is evaluated randomly, according to the cumulative
        # frequency of utilization for that appliance)
        count = 0
        maxtries = 10

        duration_indices = np.size(load_profile[load_profile > 0])
        while np.any(house_load_profile + load_profile > power_max) and count < maxtries:
        
            switch_on_index = np.size((house_load_profile + load_profile)[house_load_profile + load_profile > power_max])
            roll_step = int((duration_indices + switch_on_index)/2)

            load_profile = np.roll(load_profile, roll_step)

            count += 1
        
        # Evaluating the energy consumption from each appliance by integrating
        # the load profile over the time. Since the time is in minutes, the
        # result is divided by 60 in order to obtain Wh.
        energy[app_ID] = np.trapz(load_profile,time_sim)/60
        
        # Injecting the power demand from each appliance into the load profile of the household
        house_load_profile[:] = house_load_profile[:] + load_profile
        

    # In case the which loop failed, the instantaneous power is saturated to the
    # maximum power anyway. This may happen because the last appliance to be considered
    # is the lighting, which is of "continuous" type, therefore its profile does 
    # not depend on a cumulative frequency (they are always on) 
    house_load_profile[house_load_profile > power_max] = power_max
    
    return(house_load_profile,energy)


# # Uncomment the following lines to test the function
# import matplotlib.pyplot as plt
# from tictoc import tic, toc

# apps_availability = np.ones(17)
# day = 'wd'
# season = 's'
# dt = 1
# time = 1440
# time_sim = np.arange(0,time,dt)


# apps, apps_ID, apps_attr = datareader.read_appliances('eltdome_report.csv',';','Input')
# ec_yearly_energy, ec_levels_dict = datareader.read_enclasses('classenerg_report.csv',';','Input')
# coeff_matrix, seasons_dict = datareader.read_enclasses('coeff_matrix.csv',';','Input')

# apps_avg_lps = {}
# apps_dcs = {}
# for app in apps_ID:

#     # app_nickname is a 2 or 3 characters string identifying the appliance
#     app_nickname = apps_ID[app][apps_attr['nickname']] 

#     # app_type depends from the work cycle for the appliance: 'continuous'|'no_duty_cycle'|'duty_cycle'|
#     app_type = apps_ID[app][apps_attr['type']]

#     # app_wbe (weekly behavior), different usage of the appliance in each type of days: 'wde'|'we','wd'
#     app_wbe = apps_ID[app][apps_attr['week_behaviour']] 

#     # app_sbe (seasonal behavior), different usage of the appliance in each season: 'sawp'|'s','w','ap'
#     app_sbe = apps_ID[app][apps_attr['season_behaviour']] 

#     # Building the name of the file to be opened and read
#     fname_nickname = app_nickname
#     fname_type = 'avg_loadprof'

#     apps_avg_lps[app] = {}

#     for season in app_sbe:
#         fname_season = season

#         for day in app_wbe:
#             fname_day = day

#             filename = '{}_{}_{}_{}.csv'.format(fname_type, fname_nickname, fname_day, fname_season)
            
#             # Reading the time and power vectors for the load profile
#             data_lp = datareader.read_general(filename,';','Input')

#             # Time is stored in hours and converted to minutes
#             time_lp = data_lp[:, 0] 
#             time_lp = time_lp*60 

#             # Power is already stored in Watts, it corresponds to the load profile
#             power_lp = data_lp[:, 1] 
#             load_profile = power_lp

#             # Interpolating the load profile if it has a different time-resolution
#             if (time_lp[-1] - time_lp[0])/(np.size(time_lp) - 1) != dt: 
#                     load_profile = np.interp(time_sim, time_lp, power_lp, period = time)

#             apps_avg_lps[app][(season, day)] = load_profile

        
    # if app_type == 'duty_cycle':
    #     fname_type = 'dutycycle'
    #     filename = '{}_{}.csv'.format(fname_type, fname_nickname)
        
    #     # Reading the time and power vectors for the duty cycle 
    #     data_dc = datareader.read_general(filename, ';', 'Input')

    #     # Time is already stored in  minutes
    #     time_dc = data_dc[:, 0] 

    #     # Power is already stored in Watts, it corresponds to the duty cycle
    #     power_dc = data_dc[:, 1] 
    #     duty_cycle = power_dc
        
    #     # Interpolating the duty-cycle, if it has a different time resolution
    #     if (time_dc[-1] - time_dc[0])/(np.size(time_dc) - 1) != dt:
    #             time_dc = np.arange(time_dc[0], time_dc[-1] + dt, dt)
    #             duty_cycle = np.interp(time_dc, power_dc)

    #     apps_dcs[app] = {'time_dc': time_dc,
    #                     'duty_cycle': duty_cycle}


# appliances_data = {
#     'apps': apps,
#     'apps_ID': apps_ID,
#     'apps_attr': apps_attr,
#     'ec_yearly_energy': ec_yearly_energy,
#     'ec_levels_dict': ec_levels_dict,
#     'coeff_matrix': coeff_matrix,
#     'seasons_dict': seasons_dict,
#     'apps_avg_lps': apps_avg_lps,
#     'apps_dcs': apps_dcs,
#     }

# params = {
#     'power_max': 3000,
#     'en_class': 'D',
#     'toll': 15,
#     'devsta': 2,
#     'ftg_avg': 100
#     }
 
# house_load_profile, energy = house_load_profiler(apps_availability, day, season, appliances_data, **params)
# plt.bar(time_sim,house_load_profile,width=dt,align='edge')

# plt.show()