import numpy as np
import matplotlib.pyplot as plt
import csv
# from pulp import *
from scipy.interpolate import interp1d
from tabulate import tabulate
from pathlib import Path

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


colors = [(230, 25, 75),
        (60, 180, 75),
        (255, 225, 25),
        (0, 130, 200),
        (245, 130, 48),
        (145, 30, 180),
        (70, 240, 240),
        (240, 50, 230),
        (210, 245, 60),
        (250, 190, 212),
        (0, 128, 128),
        (220, 190, 255),
        (170, 110, 40),
        (255, 250, 200),
        (128, 0, 0),
        (170, 255, 195),
        (128, 128, 0),
        (255, 215, 180),
        (0, 0, 128),
        (128, 128, 128)]

# Transforming into rgb triplets
colors_rgb = []
for color in colors:
    color_rgb = []
    for value in color:
        color_rgb.append(float(value)/255) 
    colors_rgb.append(tuple(color_rgb))



###########################################################################################################################################################################

###############################################################################

# The basepath of the file is stored in a variable 
basepath = Path(__file__).parent

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

# Contractual power for each household (kW)
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
grid_power_max = power_max*n_hh


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


## Consumption from the aggregate of households 

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

# # Uncomment to check on energy consumption before and after interpolation between seasons and months
# yearly_consumption = 0
# for season in seasons:
#     ss = seasons[season][0]
#     for day in days:
#         dd = days[day][0]
#         n_days_day_type = days_distr[season][day]
#         yearly_consumption += np.sum(consumption_seasons[ss, :, dd])*dt*n_days_day_type

# print('Energy consumption before interpolation: {} kWh/year'.format(yearly_consumption))

# yearly_consumption = 0
# for month in months:
#     mm = months[month]['id'][0]
#     for day in days:
#         dd = days[day][0]
#         n_days_day_type = months[month]['days_distr'][day]
#         yearly_consumption += np.sum(consumption_month_day[:, mm, dd])*dt*n_days_day_type

# print('Energy consumption after interpolation: {} kWh/year'.format(yearly_consumption))



### Energy shared from the aggregate of household during one year, for all possible configurations

# Total number of configurations to be analysed
n_pv_sizes = len(pv_size_range)
n_battery_sizes = len(battery_size_range)
n_configurations = n_pv_sizes*n_battery_sizes

fixed_analysis_flag = 0
if n_configurations == 1: fixed_analysis_flag = 1

# Initializing a list where to store some relevant results to be printed as a tabulate
tab_results = []

# Initializing arrays where to store some relevant results to be plotted
if n_configurations  != 1:

    # Index of self-sufficiency in one year for each configuration
    iss_configurations = np.zeros((n_pv_sizes, n_battery_sizes))

    # Index of self-consumption in one year for each configuration
    isc_configurations = np.zeros((n_pv_sizes, n_battery_sizes))

    # Shared energy in one year for each configuration
    esh_configurations = np.zeros((n_pv_sizes, n_battery_sizes))



# A message is printed to let the user know how the simulation is proceeding 8every 25% of progress)
perc_configurations = 25
count_configurations = 0

# The user is informed that the evaluation of the shared energy is starting
message = '\nEvaluation of the energy shared by the aggregate of households for {:d} configuration(s).'.format(n_configurations)
print(message)

# Creating the input power dictionary to be passed to the method  that evaluates the shared energy
# for the configuration

input_powers_dict = {
    'pv_production_unit': pv_production_unit,
    'consumption_month_day': consumption_month_day,
    }

# Running through the different sizes for the PV

tic()
for pv_size in pv_size_range:

    # Storing the index of the PV size in the pv_size_range list (needed for properly storing the results)
    pv_index = pv_size_range.index(pv_size)

    # Storing the PV size in the technologies_dict that is passed to the methods
    technologies_dict['pv_size'] = pv_size

    # Running through the different sizes for the battery
    for battery_size in battery_size_range:

        # Storing the index of the battery size in the battery_size_range list (needed for properly storing the results)
        battery_index = battery_size_range.index(battery_size)

        # Storing the battery size in the technologies_dict that is passed to the methods
        technologies_dict['battery_size'] = battery_size

        
        ## Calling the method that computes the energy shared during one year

        results = shared_energy_evaluator(time_dict, input_powers_dict, technologies_dict, auxiliary_dict, fixed_analysis_flag)

        # Storing the results

        # Status of the optimization during the current typical day
        optimization_status = results['optimization_status']

        # Energy produced by the PV during each month (kWh/month)
        pv_production_energy = results['pv_production_energy']

        # Energy consumed by the aggregate of households during each month (kWh/month)
        consumption_energy = results['consumption_energy']

        # Excess energy fed into the grid during each month (kWh/month)
        grid_feed_energy = results['grid_feed_energy']

        # Deficit of energy purchased from the grid during each month (kWh/month)
        grid_purchase_energy = results['grid_purchase_energy']

        # Excess energy used to charge the battery during each month (kWh/month)
        battery_charge_energy = results['battery_charge_energy']

        # Deficit of energy taken from the battery during each month (kWh/month)
        battery_discharge_energy = results['battery_discharge_energy']

        # Energy shared by the aggregate of households, according to the definition provided by the Decree-Law 162/2019 (kWh/month)
        shared_energy = results['shared_energy']

        # Processing the results

        # Index of self-sufficiency in each month (%)
        index_self_suff_month = shared_energy/(consumption_energy + battery_charge_energy)*100

        # Index of self-consumption in each month (%)
        index_self_cons_month = shared_energy/(pv_production_energy + battery_discharge_energy)*100

        # Index of self-sufficiency in a year (%)
        index_self_suff_year =  np.sum(shared_energy)/np.sum(consumption_energy + battery_charge_energy)*100

        # Index of self-consumption in a year (%)
        index_self_cons_year = np.sum(shared_energy)/np.sum(pv_production_energy + battery_discharge_energy)*100


        ## Post-processing of the results 
        
        # In case of fixed size simulation for both the PV and the battery detailed results are provided
        if n_configurations == 1:

            # Storing the results in a tabular form 

            headers = ['Month', 'Week-day', 'Weekend-day', 'ISS \n(%)', 'ISC \n(%)', \
                'Shared energy \n(kWh)', 'PV production \n(kWh)', 'Consumption \n(kWh)', \
                'Grid feed \n(kWh)', 'Grid purchase \n(kWh)', 'Battery charge \n(kWh)', 'Battery discharge \n(kWh)']
            
            # Monthy results 
            for month in months:
                mm = months[month]['id'][0]
                row = [month.capitalize()]

                # Optimization status
                for day in days:
                    row = row + [optimization_status[month][day]]

                # Energies
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

                tab_results.append(row)

            # Yearly results
            tab_results.append([]*len(headers))
            tab_results.append(['Year', '/', '/', '{0:.2f}'.format(index_self_suff_year), \
                        '{0:.2f}'.format(index_self_cons_year), \
                        '{0:.1f}'.format(np.sum(shared_energy)), \
                        '{0:.1f}'.format(np.sum(pv_production_energy)), \
                        '{0:.1f}'.format(np.sum(consumption_energy)), \
                        '{0:.1f}'.format(np.sum(grid_feed_energy)), \
                        '{0:.1f}'.format(np.sum(grid_purchase_energy)), \
                        '{0:.1f}'.format(np.sum(battery_charge_energy)), \
                        '{0:.1f}'.format(np.sum(battery_discharge_energy)), \
                        ])

        # In case at least one technology is subject to a parametric analysis, only yearly quantities are shown
        else:

            # Storing the yearly values of the self-sufficiency and self-consumption indices and of the shared energy 
            # for the current configuration

            iss_configurations[pv_index, battery_index] = index_self_suff_year
            isc_configurations[pv_index, battery_index] = index_self_cons_year
            esh_configurations[pv_index, battery_index] = np.sum(shared_energy)

            # Stroing the results in a tabularr form

            # List of days in which the optimization did not work or was infeasible
            opt_did_not_work_list = []
            # opt_unnecessary_list = []
            opt_infeasible_list = []

            # Optimization status
            for month in months:
                month_nickname = months[month]['id'][1]

                for day in days:
                    day_nickname = days[day][1]

                    opt_status = optimization_status[month][day].lower().strip('\'",._- ').replace('.','').replace(' ','_')
                    if opt_status == 'opt_did_not_work': opt_did_not_work_list.append('{}, {}'.format(month_nickname.capitalize(), day_nickname))
                    # elif opt_status == 'opt_unnecessary': opt_unnecessary_list.append('{}, {}'.format(month_nickname.capitalize(), day_nickname))
                    elif opt_status == 'infeasible': opt_infeasible_list.append('{}, {}'.format(month_nickname.capitalize(), day_nickname))

            # Energies
            row = [pv_size, battery_size, '\n'.join(opt_did_not_work_list), '\n'.join(opt_infeasible_list), \
                '{0:.2f}'.format(index_self_suff_year), \
                '{0:.2f}'.format(index_self_cons_year), \
                '{0:.1f}'.format(np.sum(shared_energy)), \
                '{0:.1f}'.format(np.sum(pv_production_energy)), \
                '{0:.1f}'.format(np.sum(consumption_energy)), \
                '{0:.1f}'.format(np.sum(grid_feed_energy)), \
                '{0:.1f}'.format(np.sum(grid_purchase_energy)), \
                '{0:.1f}'.format(np.sum(battery_charge_energy)), \
                '{0:.1f}'.format(np.sum(battery_discharge_energy)), \
                ]

            tab_results.append(row)

        # The number of configurations evaluated is updated and, in case, the progress is printed
        count_configurations += 1
        if int(count_configurations/n_configurations*100) >= perc_configurations:
            print('{:d} % completed'.format(int(count_configurations/n_configurations*100)))
            perc_configurations += 25

# Header to be used in the tabulate in case of parametric analysis at least for one of the two technologies
if n_configurations != 1:
    headers = ['PV size \n(kW)', 'Battery size \n(kWh)', 'Opt. did not work in', 'Opt. infeasible in', \
                'ISS \n(%)', 'ISC \n(%)', 'Shared energy \n(kWh/year)', 'PV production \n(kWh/year)', 'Consumption \n(kWh/year)', \
                'Grid feed \n(kWh)', 'Grid purchase \n(kWh)', 'Battery charge \n(kWh)', 'Battery discharge \n(kWh)']

  
print('\n{0:d} configuration(s) evaluated in {1:.3f} s.'.format(n_configurations, toc()))

message = '\nOptimization status and results (self-sufficiency and self-consumption indices, shared energy, etc.)\n'
print(message) 
print(tabulate(tab_results, headers = headers))

    
dirname = 'Output'

# Creating an /Output/Figures folder, if not already existing
subdirname = 'Figures'

try: Path.mkdir(basepath / dirname / subdirname)
except Exception: pass

# Creating a subfolder, if not already existing
subsubdirname = '{}_{}_{}_simulation_results'.format(location, en_class, n_hh)

try: Path.mkdir(basepath / dirname / subdirname / subsubdirname)
except Exception: pass

if n_configurations != 1:

    main_size_range = pv_size_range
    lead_size_range = battery_size_range

    n_sizes_main = len(main_size_range)
    n_sizes_lead = len(lead_size_range)

    plot_specs = {
        0: {'type': 'plot', 'yaxis': 'right', 'label': 'ISS'},
        1: {'type': 'plot', 'yaxis': 'right', 'label': 'ISC'},
        2: {'type': 'bar', 'yaxis': 'left', 'label': 'Shared energy'},
        }

    fig_specs = {
        'suptitle': 'Parametric analysis',
        'xaxis_label': 'PV size (kW)',
        'yaxis_right_label': 'Performance index (%)',
        'yaxis_left_label': 'Energy (kWh/year)',
        'lead_size_name': 'Battery',
        'lead_size_uom': 'kWh',
        'yaxis_right_ylim': [0, 1.1*np.max(np.maximum(iss_configurations, isc_configurations))],
        'yaxis_left_ylim': [0, 1.1*np.max(esh_configurations)],
        }

    data = data = np.stack((iss_configurations, isc_configurations, esh_configurations), axis = 2)

    if n_pv_sizes == 1:
        main_size_range = battery_size_range
        lead_size_range = pv_size_range
        fig_specs['xaxis_label'] = 'Battery size (kWh)'
        fig_specs['lead_size_name'] = 'PV'
        fig_specs['lead_size_uom'] = 'kW'
        data = np.transpose(data, axes = (1, 0, 2))

    fig = plot.parametric_analysis(main_size_range, lead_size_range, data, plot_specs, fig_specs, **{'orientation': 'vertical'})
    
    # filename = '{}_{}_{}_{}_parametric_analysis.png'.format(location, en_class, season, n_hh)
    filename = 'parametric_analysis.png'
    fpath = basepath / dirname / subdirname / subsubdirname
            
    fig.savefig(fpath / filename) 

    plt.show()
    plt.close(fig)


if n_configurations == 1:

    pv_production_month_day = results['pv_production_month_day']
    consumption_month_day = results['consumption_month_day']
    grid_feed_month_day = results['grid_feed_month_day']
    grid_purchase_month_day = results['grid_purchase_month_day']
    battery_charge_month_day = results['battery_charge_month_day']
    battery_discharge_month_day = results['battery_discharge_month_day']
    battery_energy_month_day = results['battery_energy_month_day']
    shared_power_month_day = results['shared_power_month_day']

    # Plotting the quantities for each month

    plot_specs = {
        0: {'plot_type': 'plot', 'plot_yaxis': 'left', 'plot_label': 'pv_production', 'plot_color': colors_rgb[0], 'plot_linestyle': '-', 'plot_marker': 's'},
        1: {'plot_type': 'plot', 'plot_yaxis': 'left', 'plot_label': 'consumption', 'plot_color': colors_rgb[3], 'plot_linestyle': '-', 'plot_marker': 'o'},
        2: {'plot_type': 'plot', 'plot_yaxis': 'left', 'plot_label': 'grid_feed', 'plot_color': colors_rgb[2], 'plot_linestyle': '-', 'plot_marker': ''},
        3: {'plot_type': 'plot', 'plot_yaxis': 'left', 'plot_label': 'grid_purchase', 'plot_color': colors_rgb[4], 'plot_linestyle': '-', 'plot_marker': 's'},
        4: {'plot_type': 'plot', 'plot_yaxis': 'left', 'plot_label': 'battery_charge', 'plot_color': colors_rgb[6], 'plot_linestyle': '--'},
        5: {'plot_type': 'plot', 'plot_yaxis': 'left', 'plot_label': 'battery_discharge', 'plot_color': colors_rgb[7], 'plot_linestyle': '--'},
        6: {'plot_type': 'bar', 'plot_yaxis': 'left', 'plot_label': 'shared_power', 'plot_color': colors_rgb[9], 'plot_alpha': 0.5},
        7: {'plot_type': 'plot', 'plot_yaxis': 'right', 'plot_label': 'battery SOC', 'plot_color': colors_rgb[11], 'plot_linestyle': '--'},
    }

    fig_specs = {
        'suptitle': '\n'.join(('\nPower fluxes during two typical days', \
                    'for {} households with {} energetic class in the {} of Italy'\
                    .format(n_hh, en_class, location.capitalize()))),
        'xaxis_label': 'Time (h)',
        'yaxis_left_label': 'Power (kW)',
        'yaxis_right_label': 'SOC (%)',
    }

    for month in months:

        mm = months[month]['id'][0]

        fig_specs['title'] = month
        
        powers = np.stack((pv_production_month_day[:, mm, :],
                        consumption_month_day[:, mm, :],
                        grid_feed_month_day[:, mm, :],
                        grid_purchase_month_day[:, mm, :],
                        battery_charge_month_day[:, mm, :],
                        battery_discharge_month_day[:, mm, :],
                        shared_power_month_day[:, mm, :],
                        battery_energy_month_day[:, mm, :]/battery_size*100),
                        axis = 0)

        fig = plot.daily_profiles(time_sim, powers, plot_specs, fig_specs, **params)

        # filename = '{}_{}_{}_{}_power_fluxes.png'.format(location, en_class, month, n_hh)
        filename = '{}_power_fluxes.png'.format(month)
        fpath = basepath / dirname / subdirname / subsubdirname
        
        fig.savefig(fpath / filename)
        plt.close(fig)

plt.show()



    
    











