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
from shared_energy_evaluator import shared_energy_evaluator
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

message = '\nSimulation setup (battery)'
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

auxiliary_dict = {
    'seasons': seasons,
    'n_seasons': n_seasons,
    'months': months,
    'n_months': n_months,
    'days': days,
    'n_days': n_days,
    'days_distr': days_distr,
    }


## Time discretization
# Time is discretized in steps of one hour (according to the definition of shared energy in the Decree Law 162/2019 - "Mille proroghe")

# Total time of simulation (h) - for each typical day
time = 24

# Timestep for the simulation (h)
dt = dt_aggr/60

# Vector of time, from 00:00 to 23:59, i.e. 24 h
time_sim = np.arange(0, time, dt)
time_length = np.size(time_sim)

time_dict = {
    'time': time,
    'dt': dt,
    'time_sim': time_sim,
    'time_length': time_length,
    }



### Input data

# Maximum power from the grid (total for the aggregate of households) (kW)
grid_power_max = power_max*n_hh/1000


## Battery specification
# The battery specifications are stored in a file that can be read using the method read_param
# from the module datareader.py, that will return a dictionary

battery_specs = datareader.read_param('battery_specs.csv', ';', 'Input')

technologies_dict = {
    'grid_power_max': grid_power_max,
    'battery_specs': battery_specs,
    }


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

# fig, ax = plt.subplots(figsize = figsize)
# for month in months:
#     mm = months[month]['id'][0]
#     ax.plot(time_sim, pv_production_unit[:, mm], label = month.capitalize())

# ax.legend()


## Consumption fro the aggregate of households 

# Asking the user if detailed graphs and information about the load profile are needed to be stored as files and
# graphs (it will take some seconds more)

message = '\nEvaluation of the load profiles for the aggregate of households.'
print(message)

message = '\nDo you want to store files with detailed information about the load profiles generated and the related energy consumption?\
           \nPress \'enter\' to skip or \
           \nEnter \'ok\' to store: '
file_store_flag = input(message).strip('\'",._- ').lower()

if file_store_flag == '': file_store_flag = 0
else: file_store_flag = 1

message = '\nDo you want detailed information about the load profiles generated and the related energy consumption to be plotted and stored as figures?\
           \nPress \'enter\' to skip or \
           \nEnter \'ok\' to plot: '
fig_store_flag = input(message).strip('\'",._- ').lower()

if fig_store_flag == '': fig_store_flag = 0
else: fig_store_flag = 1

# The consumption from the aggregate of household is evaluated using the method aggregate_load_profiler
# from the module aggregate_load_profiler.py that returns a 3d-array where the load profiles for typical
# days (with a time-step of dt_aggr) are showed. The typical days are divided by season (axis = 0)  and
# day type (weekday or weekend day, axis =2). The power is given in W, therefore it has to be converted into kW
# Nota bene: only the consumption from households is evaluated (i.e. no consumption from shared commodities)

tic()  
consumption_seasons = agr_hlp(params, file_store_flag, fig_store_flag)/1000
print('\nLoad profiles evaluated in {0:.3f} s.'.format(toc()))

# The method aggregate_load_profiler computes the load profiles for eight typical days (two for each season)
# Here twelve months are considered, therefore the "seasonal" load profiles are interpolated into the months.
# To do this the represententative load profile of each season is assigned to the first month of the season
# (e.g. winter -> january). The interpolation is linear, therefore interp1d could have been used, but it does
# not allow for periodic interpolation. Therefore a for-loop in time is used, interpolating the power
# during each timestep

# Initializing a a 3d-array where where to store the profiles interpolated for each month
consumption_month_day = np.zeros((time_length, n_months, n_days))

for day in days:

    # fig, ax = plt.subplots(figsize = figsize)

    dd = days[day][0]

    for timestep in range(time_length):
        consumption_month_day[timestep, :, dd] = \
        np.interp(months_slice, seasons_slice, consumption_seasons[:, timestep, dd], period = n_months)

    # for month in months:
    #     mm = months[month]['id'][0]
    #     ax.plot(time_sim, consumption_month_day[:, mm, dd], label = month.capitalize())

    # ax.set_title(day)
    # ax.legend()

# plt.show()

input_powers_dict = {
    'pv_production_unit': pv_production_unit,
    'consumption_month_day': consumption_month_day,
    }


tab_sizes = []


# Running through the different sizes for the PV
for pv_size in pv_size_range:

    # Storing the PV size in the technologies_dict that is passed to the methods
    technologies_dict['pv_size'] = pv_size

    # Running through the different sizes for the battery
    for battery_size in battery_size_range:

        # Storing the battery size in the technologies_dict that is passed to the methods
        technologies_dict['battery_size'] = battery_size

        
        results = shared_energy_evaluator(time_dict, input_powers_dict, technologies_dict, auxiliary_dict)


        optimization_status = results['optimization_status']

        pv_production_energy = results['pv_production_energy']
        consumption_energy = results['consumption_energy']
        grid_feed_energy = results['grid_feed_energy']
        grid_purchase_energy = results['grid_purchase_energy']
        battery_charge_energy = results['battery_charge_energy']
        battery_discharge_energy = results['battery_discharge_energy']
        shared_energy = results['shared_energy']

        index_self_suff_month = shared_energy/(consumption_energy + battery_charge_energy)*100
        index_self_cons_month = shared_energy/(pv_production_energy + battery_discharge_energy)*100

        index_self_suff_year =  np.sum(shared_energy)/np.sum(consumption_energy + battery_charge_energy)*100
        index_self_cons_year = np.sum(shared_energy)/np.sum(pv_production_energy + battery_discharge_energy)*100

        # Post-processing of the results in case of fixed size simulation for both the PV and the battery 
        # (detailed results are provided)

        if len(pv_size_range) == 1 and len(battery_size_range) == 1:


        
            # Optimization status
            tab = []
            headers = ['Month', 'Week-day', 'Weekend-day', 'ISS (%)', 'ISC (%)', \
                'Shared energy (kWh)', 'PV production (kWh)', 'Consumption (kWh)', \
                'Grid feed (kWh)', 'Grid purchase (kWh)', 'Battery charge (kWh)', 'Battery discharge (kWh)']
            
            for month in months:

                mm = months[month]['id'][0]

                row = [month.capitalize()]

                for day in days:
                    row = row + [optimization_status[month][day]]

                row = row + (['{0:.2f}'.format(index_self_suff_month[mm]), \
                            '{0:.2f}'.format(index_self_cons_month[mm]), \
                            '{0:.1f}'.format(shared_energy[mm]), \
                            '{0:.1f}'.format(pv_production_energy[mm]), \
                            '{0:.1f}'.format(consumption_energy[mm]), \
                            '{0:.1f}'.format(grid_feed_energy[mm]), \
                            '{0:.1f}'.format(grid_purchase_energy[mm]), \
                            '{0:.1f}'.format(battery_charge_energy[mm]), \
                            '{0:.1f}'.format(battery_discharge_energy[mm]), \
                            ])

                tab.append(row)


            tab.append([]*len(headers))
            tab.append(['Year', '/', '/', '{0:.2f}'.format(index_self_suff_year), \
                        '{0:.2f}'.format(index_self_cons_year), \
                        '{0:.1f}'.format(np.sum(shared_energy)),
                        '', \
                        '', \
                        '', \
                        '', \
                        '', \
                        '', \
                        ])

            message = '\nOptimization status and results (self-sufficiency and self-consumption indices, shared energy)'
            print(message) 
            print(tabulate(tab, headers = headers))


        
        else:

            opt_did_not_work_list = []
            opt_unnecessary_list = []
            opt_infeasible_list = []

            for month in months:

                month_nickname = months[month]['id'][1]
                for day in days:

                    day_nickname = days[day][1]

                    opt_status = optimization_status[month][day].lower().strip('\'",._- ').replace('.','').replace(' ','_')
                    if opt_status == 'opt_did_not_work': opt_did_not_work_list.append('{}, {}'.format(month_nickname.capitalize(), day_nickname))
                    # elif opt_status == 'opt_unnecessary': opt_unnecessary_list.append('{}, {}'.format(month_nickname.capitalize(), day_nickname))
                    elif opt_status == 'infeasible': opt_infeasible_list.append('{}, {}'.format(month_nickname.capitalize(), day_nickname))

            row = [pv_size, battery_size, '\n'.join(opt_did_not_work_list), '\n'.join(opt_infeasible_list), \
                '{0:.2f}'.format(index_self_suff_year), '{0:.2f}'.format(index_self_cons_year), '{0:.1f}'.format(np.sum(shared_energy)), \
                ]

            

            tab_sizes.append(row)


headers = ['PV size (kW)', 'Battery size (kWh)', 'Opt. did not work in', 'Opt. infeasible in', \
            'ISS (%)', 'ISC (%)', 'Shared energy (kWh/year)']
message = '\nOptimization status and results (self-sufficiency and self-consumption indices, shared energy)'
print(message) 
print(tabulate(tab_sizes, headers = headers))



