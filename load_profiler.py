# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 19:46:39 2020

@author: giamm
"""

import numpy as np
import math
from scipy.interpolate import interp1d
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt

import datareader #routine created to properly read the files needed in the following 
from cumulative_frequency import cum_freq #routine created to evaluate the cumulative frequency for a given appliance, over 24h
from profile_interpolation import interp_profile #routine created to interpolate a profile

###############################################################################   
 
   
# This file contains a method that, for a given appliance (considering its 
# frequency density, its cumulative frequency, its duty cycle and its nominal
# yearly energy consumption) returns the load profile for the appliance during 
# one day (1440 min), with a resolution of 1 minute.


###############################################################################
 
def load_profiler(app, day, season, appliances_data, **params):
    
    ''' The method returns a load profile for a given appliance in a total simulation time of 1440 min, with a timestep of 1 min.
    
    Inputs:
        app - str, name of the appliance(str)
        day - str, type of day (weekday: 'wd'| weekend: 'we')
        season - str, season (summer: 's', winter: 'w', autumn or spring: 'ap')
        appliances_data - dict, various input data related to the appliances
        params - dict, simulation parameters
    
    Output:
        load_profile - 1d-array, load profile for the appliance (W) 
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
    en_class = params['en_class']

    # Tolerance and standard deviation on the appliances' duration
    toll = params['toll']
    devsta = params['devsta']

    # Average square footage of the household
    ftg_avg = params['ftg_avg']


    ## Input data for the appliances
    # Appliances' attributes, energy consumptions and user's coefficients 

    # apps is a 2d-array in which, for each appliance (rows) and attribute value is given (columns)
    apps_ID = appliances_data['apps_ID']
    
    # apps_ID is a dictionary in which, for each appliance (key), its ID number, type, weekly and seasonal behavior, class are given (value)
    apps = appliances_data['apps']

    # apps_attr is a dictionary in which the name of each attribute (value) is linked to its columns number in apps (key)
    apps_attr = appliances_data['apps_attr']
    
    # ec_yearly_energy is a 2d-array in which for each appliance, its yearly energy consumption is given for each energetic class
    ec_yearly_energy = appliances_data['ec_yearly_energy']

    # ec_levels_dict is a dictionary that links each energetic level (value) to its columns number in ec_yearly_energy
    ec_levels_dict = appliances_data['ec_levels_dict']

    # coeff_matrix is a 2d-array in which for each appliance, its coefficient k, related to user's behaviour in different seasons, is given
    coeff_matrix = appliances_data['coeff_matrix']

    # seasons_dict is a dictionary that links each season (value) to its columns number in coeff_matrix
    seasons_dict = appliances_data['seasons_dict']

    # apps_avg_lps is a dictionary where the average load profiles are stored for each appliance, according to the different seasonal and weeky behaviour
    apps_avg_lps = appliances_data['apps_avg_lps']

    # apps_dcs is a dictionary where the typical duty_cycle for "duty cycle" type appliances is stored
    apps_dcs = appliances_data['apps_dcs']
    

    ## Appliance attributes
    #  Extracting the correct data for the appliance from the applinces's attributes

    # The ID number of the appliance is stored in a variable since it will be used man times
    app_ID = apps_ID[app][apps_attr['id_number']]
    
    # T_on is the duration of the time period in which an appliance is switched on during a day (min)
    T_on = apps[app_ID, apps_attr['time_on'] - (len(apps_attr) - np.size(apps, 1))]
    
    # energy is the yearly energy consumption, according to the energy label (kWh/year)
    energy = ec_yearly_energy[app_ID, ec_levels_dict[en_class]]
    
    # kk is the coefficient accounting for the user's behaviour in different season, for each appliance(-)
    kk = coeff_matrix[app_ID, seasons_dict[season]]

    # Return a vector of zeros if the appliance is not used at all during the season
    if kk == 0: return(np.zeros((np.shape(time_sim))))

    # k_ftg is a coefficient that adapts the consumption from lux to the actual footage of the household
    k_ftg = ftg_avg/100 #(m2/m2)

    ## Nominal power consumption     
    # The nominal power is evaluated, according to the yearly energy 
    # consumption (kWh/year) and the time-period in which the appliance is 
    # switched on during a day (min/day). The first one is given in kWh/year
    # and it is multiplied by 1000 and divided by 365 in order to get Wh/day.
    # The second one is given in minutes and it is divided by 60 to get hours.
    
    # An if statement is used to avoid negative power when T_on = -1
    if T_on == -1:
        power = 0
    else:
        power = (energy*1000/365)/(T_on/60) #(W)
    
    ## Input data for the appliance
    # Selecting the correct key where to find the average daily load profile in the dictionary
    
    # app_nickname is a 2 or 3 characters string identifying the appliance
    app_nickname = apps_ID[app][apps_attr['nickname']] 

    # app_type depends from the work cycle for the appliance: 'continuous'|'no_duty_cycle'|'duty_cycle'|
    app_type = apps_ID[app][apps_attr['type']]

    # app_wbe (weekly behavior), different usage of the appliance in each type of days: 'wde'|'we','wd'
    app_wbe = apps_ID[app][apps_attr['week_behaviour']] 

    # app_sbe (seasonal behavior), different usage of the appliance in each season: 'sawp'|'s','w','ap'
    app_sbe = apps_ID[app][apps_attr['season_behaviour']] 

    # Default choice (no different behaviour for different seasons): 
    # if the appliance has got different profiles in different seasons, this will be changed
    key_season = 'sawp' 
    if len(app_sbe) > 1: key_season = season

    # Default choice (no different behaviour for different types of day):
    # if the appliance has got different profiles in different days of the week, this will be changed
    key_day = 'wde' 
    if len(app_wbe) > 1: key_day = day

    avg_load_profile = apps_avg_lps[app][(key_season, key_day)]
    


    ### Different routines according to the appliance's type
    

    ## Routine to be followed for appliances which belong to "continuous" type (ac and lux)
    # The load profile has already been evaluated therefore it only needs to be loaded and properly managed
    
    if app_type == 'continuous':
        
        load_profile = avg_load_profile

        # Lighting: the load profile is taken as it is, since no information about the different 
        # yearly energy consumption are available. The value is just adjusted to the users' seasonal
        # behavior (kk) and house's footage (ftg_avg) 
        if app_nickname == 'lux':
            activate_kk = 0
            if activate_kk == 0: kk = 1  
            load_profile = load_profile*kk*k_ftg 

        # Other types of continuous appliances: the load profile is adjusted to the yearly energy consumption
        else:
            activate_kk = 1
            if activate_kk == 0: kk = 1  
            load_profile = load_profile/(np.trapz(load_profile, time_sim)/(time_sim[-1] - time_sim[0]))*power*kk


        return(load_profile)
   
     
    ## Routine to be followed for appliances which belong to "duty cycle" or "uniform" types
    # The load profile has to be evaluated according to the time-period in which they are switched on (T_on)
    # and to the frequency density of their usage during the day
     
    # Loading the duty-cycle, for those appliances which have one
    if app_type == 'duty_cycle':
        
        time_dc = apps_dcs[app]['time_dc']
        duty_cycle = apps_dcs[app]['duty_cycle']

        # Adjusting the duty cycle to the actual nominal power 
        # and the users' habits (coefficient kk, varying according to the season)
        activate_kk = 0
        if activate_kk == 0: kk = 1
        duty_cycle = duty_cycle/(np.trapz(duty_cycle, time_dc)/(time_dc[-1] - time_dc[0]))*power*kk
    
    # Building a uniform duty-cycle for those appliances which use a constant power,
    # according to the time-period in a day in which they are used (T_on)  
    if app_type == 'uniform':
       
        # The following procedure is performed only if the appliance is not used continously during the day
        if T_on != time: 
        
            # A randomly variating (in a normal distribution) duration is 
            # assumed for appliances which don't have a duty-cycle
            d = T_on
            lim_up = d+(toll*(d)/100)        
            lim_low = d-(toll*(d)/100)
            lim_up , lim_low = np.clip((lim_up, lim_low), 0, time) 
            T_on = 0
            
            while (T_on > lim_up or T_on < lim_low):
                T_on = int(np.random.normal(d, devsta))

        # Adjusting T_on to the seasonal coefficient for the user's behaviour
        T_on = T_on*kk

        # Creating the time and power vectors for the (uniform) duty-cycle      
        time_dc = np.linspace(0, T_on, int(T_on/dt + 1))
        duty_cycle = np.ones(np.shape(time_dc))
        
        # Giving the duty cycle the same shape as "duty_cycle" type appliances, 
        # meaning that the first value for the power is 0
        duty_cycle[0] = 0
     
        # Adjusting the duty cycle to the actual nominal power 
        duty_cycle = duty_cycle/(np.trapz(duty_cycle, time_dc)/(time_dc[-1] - time_dc[0]))*power
    

    ## Usage probability distributions
    
    # Selecting a time instant from the usage's frquency distribution of the appliance. The latter
    # is equal to the average daily load profile (a normalization is perfoemd since the latter is in W)
    freq_dens = apps_avg_lps[app][(key_season, key_day)]
    
    # Evaluating the cumulative frquency of appliance'usage in one day
    # cumfreq = cum_freq(time_sim, freq_dens)
    cum_freq = cumtrapz(freq_dens, time_sim, initial = 0)
    cum_freq = cum_freq/np.max(cum_freq)
 

    ## Switch-on instant
    # Selecting a random istant in which the appliances starts working (according to the cumulative frequency) 
    # and using its duty-cycle (uniform for those appliances which don't have a duty-cycle) to create the load profile
    
    # Selecting a random instant when to make the appliance start its cycle, 
    # according to the frequency density and the cumulative frequency
    random_probability = np.random.rand()
    
    # Evaluating a time instant at which the appliance starts its cycle through the cumulative frequency,
    # extracting it from a probability distribution that follows the frequency density of the appliance's usage
    switch_on_instant = time_sim[cum_freq >= random_probability][0]
    switch_on_index = int(np.where(time_sim == switch_on_instant)[0]) 
    

    ## Building the load profile using the appliance's duty-cycle.
    # This is done by initializing a vector of zeros of the same length as time_sim
    # then injecting the duty-cycle at the begininng and shifting it to the switch-on
    # instant using np.roll
    load_profile = np.zeros(np.shape(time_sim))
    load_profile[:len(duty_cycle)] = duty_cycle
    load_profile = np.roll(load_profile, switch_on_index)

    return (load_profile)


# ## Uncomment the following lines to test the function (comment def and return)
# # app = 'vacuum_cleaner'
# # app = 'air_conditioner'
# # app = 'electric_oven'
# # app = 'microwave_oven'
# # app = 'fridge'
# # app = 'freezer'
# app = 'washing_machine'
# # app = 'dish_washer'
# # app = 'tumble_drier'
# # app = 'electric_boiler'
# # app = 'hifi_stereo'
# # app = 'dvd_reader'
# # app = 'tv'
# # app = 'iron'
# # app = 'pc'
# # app = 'laptop'
# # app = 'lighting'

# day = 'we'
# season = 's'

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
#     'en_class': 'D',
#     'toll': 15,
#     'devsta': 2,
#     'ftg_avg': 100,
# }



# dt = 1 
# time = 1440 
# time_lp = np.arange(0,time,dt)
# load_profile=load_profiler(app,day,season,appliances_data, **params)
# # plt.bar(time_lp,load_profile,width=time_lp[-1]/np.size(time_lp))
# plt.plot(time_lp,load_profile)

# switch_on_index = int(np.where(load_profile > 0)[0][0]) 
# print('switch_on_time{} min - index: {}'.format(time_lp[switch_on_index], switch_on_index))
# print('power at switch-on time: {} W'.format(load_profile[switch_on_index]))

# print('Daily consumption ({}): {} Wh/day'.format(app, np.trapz(load_profile, time_lp)/60))
# print('Yearly consumption ({}): {} kWh/year'.format(app, np.trapz(load_profile, time_lp)/60*365/1000))

# plt.show()


#####################################################################################

    # # Building the name of the file to be opened and read
    # fname_nickname = app_nickname
    
    # # Default choice (no different behaviour for different types of day):
    # # if the appliance has got different profiles in different days of the week, this will be changed
    # fname_day = 'wde' 
    # if len(app_wbe) > 1: fname_day = day

    # # Default choice (no different behaviour for different seasons): 
    # # if the appliance has got different profiles in different seasons, this will be changed
    # fname_season = 'sawp' 
    # if len(app_sbe) > 1: fname_season = season

######### Old way of loading load profile files (type == 'continuous')
        # fname_type = 'avg_loadprof'
        # filename = '{}_{}_{}_{}.csv'.format(fname_type, fname_nickname, fname_day, fname_season)
        
        # # Reading the time and power vectors for the load profile
        # data_lp = datareader.read_general(filename,';','Input')

        # # Time is stored in hours and converted to minutes
        # time_lp = data_lp[:, 0] 
        # time_lp = time_lp*60 

        # # Power is already stored in Watts, it corresponds to the load profile
        # power_lp = data_lp[:, 1] 
        # load_profile = power_lp
        
        # # Interpolating the load profile if it has a different time-resolution
        # if (time_lp[-1] - time_lp[0])/(np.size(time_lp) - 1) != dt: 
        #     # load_profile = interp_profile(time_lp, power_lp, time_sim)
        #     load_profile = np.interp(time_sim, time_lp, power_lp, period = time)

######### Old way of loading duty cycle files (type == 'duty_cycle')
        # fname_type = 'dutycycle'
        # filename = '{}_{}.csv'.format(fname_type, fname_nickname)
        
        # # Reading the time and power vectors for the duty cycle 
        # data_dc = datareader.read_general(filename, ';', 'Input')

        # # Time is already stored in  minutes
        # time_dc = data_dc[:, 0] 

        # # Power is already stored in Watts, it corresponds to the duty cycle
        # power_dc = data_dc[:, 1] 
        # duty_cycle = power_dc
        
        # # Interpolating the duty-cycle, if it has a different time resolution
        # if (time_dc[-1] - time_dc[0])/(np.size(time_dc) - 1) != dt:
        #     time_dc = np.arange(time_dc[0], time_dc[-1] + dt, dt)
        #     duty_cycle = np.interp(time_dc, power_dc)
        #     # f = interp1d(time_dc, power_dc, 
        #     #     kind = 'linear', bounds_error = False, fill_value = 'extrapolate')
        #     # time_dc = np.arange(time_dc[0], time_dc[-1] + dt, dt)
        #     # duty_cycle = f(time_dc)

######### Old way of loading frequency density files (type != 'continuous')
        # # Loading the average daily load profile during one day
        # # It is used as a frequency distribution for the appliance's use
        # # and it is managed in order to get a cumulative frequency
        # fname_type = 'avg_loadprof'
        # filename = '{}_{}_{}_{}.csv'.format(fname_type, fname_nickname, fname_day, fname_season)

        # # Reading the time and power vectors for the frequency density
        # data_fq = datareader.read_general(filename, ';', 'Input')

        # # Time is stored in hours and converted into minutes
        # time_fq = data_fq[:, 0] 
        # time_fq = time_fq*60 

        # # Power is already stored in Watts, it corresponds to the frequency density
        # power_fq = data_fq[:, 1] 
        # freq_dens = power_fq
        
        # # Interpolating the frqeuency density, if it has a different time resolution
        # if (time_fq[-1])/(np.size(time_fq) - 1) != dt:
        #         # freq_dens = interp_profile(time_fq, power_fq, time_sim)
        #         freq_dens = np.interp(time_sim, time_fq, power_fq, period = time)