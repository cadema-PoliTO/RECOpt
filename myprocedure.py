import numpy as np
import matplotlib.pyplot as plt
import csv
# from pulp import *
from scipy.interpolate import interp1d
from tabulate import tabulate

import parameters_input as inp
import plot_generator as plot
import datareader
from tictoc import tic, toc
from aggregate_load_profiler import aggregate_load_profiler as agr_hlp
from battery_optimization import battery_optimization







## Plot parameters

# Default parameters
def_params = {
'time_scale': 'h',
'power_scale': 'kW',
'energy_scale': 'MWh',
'figsize': (297/25.4 , 420/25.4),
'orientation': 'horizontal',
'font_small': 14,
'font_medium': 16,
'font_large': 18,
}

# The parameters that are not specified when the function is called are set to the default value
params = {}
for param in def_params: 
    if param not in params: params[param] = def_params[param]

# Figure setup: figure size and orientation, font-sizes 
figsize = params['figsize']
orientation = params['orientation']

if orientation == 'horizontal': figsize = figsize[::-1]

fontsize_title = params['font_large']
fontsize_legend = params['font_medium']
fontsize_labels = params['font_medium']
fontsize_text = params['font_medium']
fontsize_ticks = params['font_small']
fontsize_pielabels = params['font_small']



###########################################################################################################################################################################

### Parameters and simulation setup


## Parameters 

# Parameters are enter from keyboard by the user. This operation is performed using
# the method parameters_input from parameters_input.py (inp) that returns a dictionary.

params = inp.parameters_input()

# The values of the parameters that are needed here are taken from params

# Number of households considered (-)
n_hh = params['n_hh'] 

# Geographical location: 'north' | 'centre' | 'south'
location = params['location']

# Energetic class of the appiances: 'A+++' | 'A++' | 'A+' | 'A' | 'B' | 'C' | 'D'
en_class = params['en_class']

# Contractual power for each household (W)
power_max = params['power_max']

# Time-step used to aggregated the results (min): 1 | 5 | 10 | 15 | 10 | 30 | 45 | 60
dt_aggr = params['dt_aggr']


## Type of simulation
# The type of simulation is chosen (fixed size or paramtric analysis) both for the pv
# and the battery and the ranges for the sizes are specified as well. This is done using the method
# simulation_setup from parameters_input.py (inp) that returns a dictionary and a list.

# Photovoltaic (PV)
pv_setup, pv_size_range = inp.simulation_setup('PV')

# Battery 
battery_setup, battery_size_range = inp.simulation_setup('battery')


## Showing the users all the updated parameters
   
message = '\nThe parameters for the simulation and the simulation setup are now set as follows'
print(message)

# Simulation parameters
tab = []
for param in params:
    row = [param, params[param]]
    tab.append(row)

message = '\nSimulation parameters'
print(message) 
print(tabulate(tab, headers=['Parameter', 'Value']))

# Simulation setup (PV)
tab = []
for param in pv_setup:
    row = [param, pv_setup[param]]
    tab.append(row)

message = '\nSimulation setup (PV)'
print(message) 
print(tabulate(tab, headers=['Parameter', 'Value']))

# Simulation setup (battery)
tab = []
for param in battery_setup:
    row = [param, battery_setup[param]]
    tab.append(row)

message = '\nSimulation parameters'
print(message) 
print(tabulate(tab, headers=['Parameter', 'Value']))


## Building seasons, months and week dictionaries 
# This is done in order to explore all the seasons and, for each season two types of days (weekday and weekend)    
 
seasons = {
    'winter': (0, 'w'),
    'spring': (1, 'ap'),
    'summer': (2, 's'),
    'autumn': (3, 'ap')
    }

days = {
    'week-day': (0, 'wd'),
    'weekend-day': (1, 'we')
    }

# While the routine that evaluated the load profiles works with typical profiles for each season, 
# the optimization procedure works with tyipcal profiles for each month.  
# A reference year is considered, in which the first day (01/01) is a monday. 
# Therefore, conventionally considering that winter lasts from january to march 
# spring from april to june, summer from july to september and autumn from october
# to december, each month has got the following number of weekdays and weekend days.

# Creating a dictionary for the months
months = {
    'january': {'id': (0, 'jan'), 'season': 'winter', 'days_distr': {'week-day': 23, 'weekend-day': 8}},
    'february': {'id': (1, 'feb'), 'season': 'winter', 'days_distr': {'week-day': 20, 'weekend-day': 8}},
    'march': {'id': (2, 'mar'), 'season': 'winter', 'days_distr': {'week-day': 22, 'weekend-day': 9}},
    'april': {'id': (3, 'apr'), 'season': 'spring', 'days_distr': {'week-day': 21, 'weekend-day': 9}},
    'may': {'id': (4, 'may'), 'season': 'spring', 'days_distr': {'week-day': 23, 'weekend-day': 8}},
    'june': {'id': (5, 'jun'), 'season': 'spring', 'days_distr': {'week-day': 21, 'weekend-day': 9}},
    'july': {'id': (6, 'jul'), 'season': 'summer', 'days_distr': {'week-day': 22, 'weekend-day': 9}},
    'august': {'id': (7, 'aug'), 'season': 'summer', 'days_distr': {'week-day': 23, 'weekend-day': 8}},
    'september': {'id': (8, 'sep'), 'season': 'summer', 'days_distr': {'week-day': 20, 'weekend-day': 10}},
    'october': {'id': (9, 'oct'), 'season': 'autumn', 'days_distr': {'week-day': 23, 'weekend-day': 8}},
    'november': {'id': (10, 'nov'), 'season': 'autumn', 'days_distr': {'week-day': 22, 'weekend-day': 8}},
    'december': {'id': (11, 'dec'), 'season': 'autumn', 'days_distr': {'week-day': 21, 'weekend-day': 10}},
    }

# The days distribution in the seasons can be evaluated as well
days_distr = {}
for month in months:
        season = months[month]['season']
        
        if season not in days_distr: days_distr[season] = {'week-day': months[month]['days_distr']['week-day'],
                                                           'weekend-day': months[month]['days_distr']['weekend-day']}

        else: 
                days_distr[season]['week-day'] += months[month]['days_distr']['week-day']
                days_distr[season]['weekend-day'] += months[month]['days_distr']['weekend-day']

# Storing some useful quantities
n_months = len(months)
n_seasons = len(seasons)
n_days = len(days)

# Creating two lists to properly slice the consumption data when interpolating from seasonal to monthly values
months_slice = np.arange(0, n_months)
seasons_slice = np.arange(0, n_months, int(n_months/n_seasons))


## Time discretization
# Time is discretized in steps of one hour (according to the definition of shared energy in the Decree Law 162/2019 - "Mille proroghe")

# Total time of simulation (h) - for each typical day
time = 24

# Timestep for the simulation (h)
dt = 1

# Vector of time, from 00:00 to 23:59, i.e. 24 h
time_sim = np.arange(0, time, dt)
time_length = np.size(time_sim)



### Input data

# Maximum power from the grid (total for the aggregate of households) (kW)
grid_power_max = power_max*n_hh/1000


## Battery specification
# The battery specifications are stored in a file that can be read using the method read_param
# from the module datareader.py, that will return a dictionary

battery_specs = datareader.read_param('battery_specs.csv', ';', 'Input')


## Unit production from the photovoltaic installation (kWh/h/kWp)
# The unit production during each hour (kWh/h/kWp) from the photovoltaic installation can 
# be read using the method read_general from the module datareader.py, that returns a 2d-array
# containing the time vector on the first columns and the unit production during each hour in each
# month in the other columns

data_pv = datareader.read_general('pv_production_unit.csv', ';', 'Input')

time_pv = data_pv[:, 0]
pv_production_unit = data_pv[:, 1:]

# Interpolating the unit pv production, if it has a different time resolution
if (time_pv[-1] - time_pv[0])/(np.size(time_pv) - 1) != dt:
    
    f_pv = interp1d(time_pv, pv_production_unit, kind = 'linear', axis = 0, fill_value = 'extrapolate')
    pv_production_unit = f(time_sim)

fig, ax = plt.subplots(figsize = figsize)
for month in months:
    mm = months[month]['id'][0]
    ax.plot(time_sim, pv_production_unit[:, mm], label = month.capitalize())

ax.legend()


## Consumption fro the aggregate of households
# The consumption from the aggregate of household is evaluated using the method aggregate_load_profiler
# from the module aggregate_load_profiler.py that returns a 3d-array where the load profiles for typical
# days (with a time-step of dt_aggr) are showed. The typical days are divided by season (axis = 0)  and
# day type (weekday or weekend day, axis =2). The power is given in W, therefore it has to be converted into kW
# Nota bene: only the consumption from households is evaluated (i.e. no consumption from shared commodities)      
    
consumption_seasons = agr_hlp(params)/1000

# The method aggregate_load_profiler computes the load profiles for eight typical days (two for each season)
# Here twelve months are considered, therefore the "seasonal" load profiles are interpolated into the months.
# To do this the represententative load profile of each season is assigned to the first month of the season
# (e.g. winter -> january). The interpolation is linear, therefore interp1d could have been used, but it does
# not allow for periodic interpolation. Therefore a for-loop in time is used, interpolating the power
# during each timestep

# Initializing a a 3d-array where where to store the profiles interpolated for each month
consumption_month_day = np.zeros((time_length, n_months, n_days))



for day in days:

    fig, ax = plt.subplots(figsize = figsize)

    dd = days[day][0]

    for timestep in time_sim:
        consumption_month_day[timestep, :, dd] = \
        np.interp(months_slice, seasons_slice, consumption_seasons[:, timestep, dd], period = n_months)

    for month in months:
        mm = months[month]['id'][0]
        ax.plot(time_sim, consumption_month_day[:, mm, dd], label = month.capitalize())

    ax.set_title(day)
    ax.legend()

plt.show()


 



# # ## Solving
# # # Mixed Integer Linear Problem (MILP) Optmiziation for the battery control
# pv_energy = np.zeros((n_months, n_days))
# feed_energy = np.zeros((n_months, n_days))
# purchase_energy = np.zeros((n_months, n_days))
# consumption_energy = np.zeros((n_months, n_days))
# battery_charge_energy = np.zeros((n_months, n_days))
# battery_discharge_energy = np.zeros((n_months, n_days))
# shared_energy = np.zeros((n_months, n_days))


# for month in months:

#     mm = months[month]['id'][0]

#     pv_production = np.zeros((time_length, n_days)) 
#     ue_consumption = np.zeros((time_length, n_days)) 
#     net_load = np.zeros((time_length, n_days)) 
#     pv_available = np.zeros((time_length, n_days))
#     battery_charge= np.zeros((time_length, n_days))
#     battery_discharge = np.zeros((time_length, n_days))
#     grid_feed = np.zeros((time_length, n_days))
#     grid_purchase = np.zeros((time_length, n_days))
#     battery_energy = np.zeros((time_length, n_days))
#     shared_power = np.zeros((time_length, n_days))


#     # # Creating a figure with multiple subplots, with two rows (one for each type of day)
#     # fig, ax = plt.subplots(2, 1, sharex = False, sharey = False, figsize = figsize)
    
#     # # suptitle = 'Results in {}'.format(month)
#     # # fig.suptitle(suptitle, fontsize = fontsize_title, fontweight = 'bold')
#     # fig.subplots_adjust(left = 0.1, bottom = 0.1, right = 0.9, top = 0.85, wspace = None, hspace = 0.3)


#     for day in days:

#         dd = days[day][0]

#         pv_production[:, dd]  = pv_production_month[:, mm]  
#         ue_consumption[:, dd] = ue_consumption_month_day[:, mm, dd]

#         pv_available[:, dd] = pv_production[:, dd] - ue_consumption[:, dd]
#         net_load[pv_available[:, dd] < 0, dd] = -pv_available[pv_available[:, dd] < 0, dd]
#         pv_available[pv_available[:, dd] < 0, dd]= 0
        
#         grid_feed[:, dd], grid_purchase[:, dd], battery_charge[:, dd], battery_discharge[:, dd], battery_energy[:, dd] = \
#         battery_optimization(pv_available[:, dd], net_load[:, dd], battery_capacity, battery_specs, grid_power, time_length, time_sim, dt)
    
#         shared_power = np.minimum((pv_production + battery_discharge), (ue_consumption + battery_charge))

#         number_of_days = months[month]['days_distr'][day]

#         pv_energy[mm][dd] = np.sum(pv_production)*dt*number_of_days
#         feed_energy[mm][dd] = np.sum(grid_feed)*dt*number_of_days
#         purchase_energy[mm][dd] = np.sum(grid_purchase)*dt*number_of_days
#         consumption_energy[mm][dd] = np.sum(ue_consumption)*dt*number_of_days
#         battery_charge_energy[mm][dd] = np.sum(battery_charge)*dt*number_of_days
#         battery_discharge_energy[mm][dd] = np.sum(battery_discharge)*dt*number_of_days
#         shared_energy[mm][dd] = np.sum(shared_power)*dt*number_of_days

#         tol = 1e-4
#         if (np.any(net_load[:, dd] + grid_feed[:, dd] + battery_charge[:, dd] - pv_available[:, dd] - grid_purchase[:, dd] - battery_discharge[:, dd] > tol)): print('ops, node')
#         if (np.any(grid_feed[:, dd] > grid_power + tol) or np.any(grid_feed[:, dd] > grid_power + tol)): print('ops, grid power')
#         if (np.any(grid_feed[:, dd] > pv_available[:,dd] + tol)): print('ops, grid feed')
#         if (np.any(battery_charge[:, dd] > pv_available[:,dd] + tol)): print('ops, battery_charge')


# pv_energy_month = np.sum(pv_energy, axis = 1)
# feed_energy_month = np.sum(feed_energy, axis = 1)
# purchase_energy_month = np.sum(purchase_energy, axis = 1)
# consumption_energy_month= np.sum(consumption_energy, axis = 1)
# battery_charge_energy_month = np.sum(battery_charge_energy, axis = 1)
# battery_discharge_energy_month = np.sum(battery_discharge_energy, axis = 1)
# shared_energy_month = np.sum(shared_energy, axis = 1)

# index_self_suff_month = shared_energy_month/consumption_energy_month*100
# index_self_cons_month = shared_energy_month/pv_energy_month*100

# index_self_suff_year =  np.sum(shared_energy_month)/np.sum(consumption_energy_month)*100
# index_self_cons_year = np.sum(shared_energy_month)/np.sum(pv_energy_month)*100

# [print('{0:s}: ISS: {1:.2f}%, ISC: {2:.2f}%'.format(month.capitalize(), index_self_suff_month[months[month]['id'][0]], index_self_cons_month[months[month]['id'][0]])) for month in months]

# print('\nYear: ISS: {0:.2f}%, ISC: {1:.2f}%'.format(index_self_suff_year, index_self_cons_year))

# #         ax[dd].plot(time_sim + dt/2, pv_available[:, dd], label = 'pv_avaialble')
# #         ax[dd].plot(time_sim + dt/2, net_load[:, dd], label = 'net_load')
# #         ax[dd].plot(time_sim + dt/2, grid_feed[:, dd], label = 'grid_feed')
# #         ax[dd].plot(time_sim + dt/2, grid_purchase[:, dd], label = 'grid_purchase')
# #         ax[dd].plot(time_sim + dt/2, battery_charge[:, dd], label = 'battery_charge')
# #         ax[dd].plot(time_sim + dt/2, battery_discharge[:, dd], label = 'battery_discharge')

# #         ax[dd].bar(time_sim, shared_power[:, dd], color = 'k', width = dt, align = 'edge', label = 'shared power', alpha = 0.3)
        
# #         title = '{}, {}'.format(month.capitalize(), day)
# #         ax[dd].set_title(title, fontsize = fontsize_title)

# #         axtw = ax[dd].twinx()
# #         axtw.bar(time_sim, battery_energy[:, dd]/battery_capacity*100, color = 'y', width = dt, align = 'edge', label = 'shared power', alpha = 0.3)
        
# #     ##
# #     # Making the figure look properly
# #     for axi in ax.flatten():
# #         axi.set_xlabel('Time ({})'.format('h'), fontsize = fontsize_labels)
# #         axi.set_ylabel('Power ({})'.format('kW'), fontsize = fontsize_labels)
# #         axi.set_xlim([time_sim[0], time_sim[-1]])
# #         # Set one tick each hour on the x-axis
# #         axi.set_xticks(list(time_sim[: : int(dt)]))
# #         axi.tick_params(axis ='both', labelsize = fontsize_ticks)
# #         axi.tick_params(axis ='x', labelrotation = 0)
# #         axi.grid()
# #         axi.legend(loc = 'upper left', fontsize = fontsize_legend, ncol = 2)
    

# # plt.show()


