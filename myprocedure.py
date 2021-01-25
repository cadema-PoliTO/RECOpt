import numpy as np
#from numpy import *
import matplotlib.pyplot as plt
from pulp import *
import csv
# import numpy_financial as npf
from tictoc import tic, toc
from scipy.interpolate import interp1d
from aggregate_load_profiler import aggregate_load_profiler as agr_hlp
from battery_optimization import battery_optimization
from parameters_input import parameters_input
import plot_generator as plot





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

params = parameters_input()

# Number of households considered (-)
n_hh = params['n_hh'] 

# Geographical location: 'north' | 'centre' | 'south'
location = params['location']

# Energetic class of the appiances: 'A+++' | 'A++' | 'A+' | 'A' | 'B' | 'C' | 'D'
en_class = params['en_class']

# Time-step used to aggregated the results (min): 1 | 5 | 10 | 15 | 10 | 30 | 45 | 60
dt_aggr = params['dt_aggr']

# Contractual power for each household (W)
power_max = params['power_max']




# message = 'Please, specify the range of possible values for the PV size\nEnter 0 or press enter to skip and let the programme decide the range'
# print(message)
# pv_size_min = input('Minimum size: ')
# pv_size_max = input('Maximum size: ')

# message = 'Please, specify the range of possible values for the battery size\nEnter 0 or press enter to skip and let the programme decide the range'
# print(message)
# battery_size_min = input('Minimum size: ')
# battery_size_max = input('Maximum size: ')


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



n_months = len(months)
n_seasons = len(seasons)
n_days = len(days)

months_slice = np.arange(0, n_months)
seasons_slice = np.arange(0, n_months, int(n_months/n_seasons))



grid_power = power_max*n_hh

pv_size = 10          # Potenza di picco PV (kW))
battery_size = 8       # Capacità batteria (kWh)

## Battery specifications
battery_capacity = battery_size
battery_specifications_file = np.loadtxt("Battery_spec.txt")  
battery_specs = {
    'SOC_min': battery_specifications_file[0],
    'SOC_max': battery_specifications_file[1],
    't_cd_min': battery_specifications_file[2],
    'eta_charge': battery_specifications_file[3],
    'eta_discharge': battery_specifications_file[4],
    'eta_self_discharge': battery_specifications_file[5],
    }

pv_production_unit = np.loadtxt("PPV.txt")                  # Produzione fotovoltaica (per 1 kWp)
pv_production_month = pv_production_unit*pv_size                  # Produzione fotovoltaica scalata alla taglia
    
       
ue_seasons = agr_hlp(params)                   # Carico elettrico delle utenze nella configurazione (kW)




## Time discretization
# Time is discretized in steps of one hour (according to the definition of shared energy in the Decree Law 162/2019 - "Mille proroghe")

# Total time of simulation (h) - for each typical day
time = 24

# Timestep for the simulation (h)
dt = 1

# Vector of time, from 00:00 to 23:59, i.e. 24 h
time_sim = np.arange(0, time, dt)
time_length = np.size(time_sim)

# seasons_n = 4
# months_n = 12 
# days_n = 2

# seasons = np.arange(0, months_n, int(months_n/seasons_n))    
# months = np.arange(0, months_n)
# days = np.arange(0, days_n)

# print(seasons)
# print(months)


## Consumption from the households
# The method aggregate_load_profiler computes the load profiles for eight typical days (a weekday and a weekend day
# for each season), that are aggregated with a specified resolution. Here twelve months are considered, therefore
# the load profiles for the seasons are interpolated assigning the load profile of each season to the first month

# Initializing a a 3d-array where where to store the profiles interpolated for each month
ue_consumption_month_day = np.zeros((time_length, n_months, n_days))

for day in days:

    dd = days[day][0]


    for timestep in time_sim:
        ue_consumption_month_day[timestep, :, dd] = \
        np.interp(months_slice, seasons_slice, ue_seasons[:, timestep, dd], period = n_months)/1000


# print(np.shape(ue_consumption))

    



 

# ## Powers to be considered
# # Initializing the vectors
# pv_available = np.zeros((time_length, n_months))            # Power from PV tht is not indirectly self-consumed (kW)
# ue_net_load = np.zeros((time_length, n_months))             # Net load for the users part of the configuration (kW)
# battery_charge= np.zeros((time_length, n_months))           # Power that charges the batter (kW)
# battery_discharge = np.zeros((time_length, n_months))       # Power discharged from the battery (kW)
# grid_feed = np.zeros((time_length, n_months))               # Power fed into the grid (kW)
# grid_purchase = np.zeros((time_length, n_months))           # Power purchased from the grid (kW)
# battery_energy = np.zeros((time_length, n_months))          # Energy stored in the battery (kWh)

# # Evaluating the actual values

# ue_consumption = np.transpose((ue_consumption[:, :, 0]/1000))

# pv_available = pv_production - ue_consumption
# ue_net_load[pv_available < 0] = -pv_available[pv_available < 0]
# pv_available[pv_available < 0] = 0




# ## Showing some data

# month = 7

# fig, ax = plt.subplots(figsize = figsize)

# ax.bar(time_sim, ue_consumption[:, month], color = 'b', align = 'edge', label = 'ue', alpha = 0.5)
# ax.plot(time_sim + dt/2, ue_consumption[:, month], color = 'b')
# ax.bar(time_sim, pv_production[:, month], color = 'r', align = 'edge', label = 'pv_prod', alpha = 0.5)
# ax.plot(time_sim + dt/2, pv_production[:, month], color = 'r')
# ax.bar(time_sim, pv_available[:, month], color = 'k', align = 'edge', label = 'pv_ava', alpha = 0.5)
# ax.plot(time_sim + dt/2, pv_available[:, month], color = 'k')
# ax.bar(time_sim, net_load[:, month], color = 'g', align = 'edge', label = 'net_load', alpha = 0.5)
# ax.plot(time_sim + dt/2, net_load[:, month], color = 'g')

# ax.legend()
# plt.show(fig)


# ## Solving
# # Mixed Integer Linear Problem (MILP) Optmiziation for the battery control
pv_energy = np.zeros((n_months, n_days))
feed_energy = np.zeros((n_months, n_days))
purchase_energy = np.zeros((n_months, n_days))
consumption_energy = np.zeros((n_months, n_days))
battery_charge_energy = np.zeros((n_months, n_days))
battery_discharge_energy = np.zeros((n_months, n_days))
shared_energy = np.zeros((n_months, n_days))


for month in months:

    mm = months[month]['id'][0]

    pv_production = np.zeros((time_length, n_days)) 
    ue_consumption = np.zeros((time_length, n_days)) 
    net_load = np.zeros((time_length, n_days)) 
    pv_available = np.zeros((time_length, n_days))
    battery_charge= np.zeros((time_length, n_days))
    battery_discharge = np.zeros((time_length, n_days))
    grid_feed = np.zeros((time_length, n_days))
    grid_purchase = np.zeros((time_length, n_days))
    battery_energy = np.zeros((time_length, n_days))
    shared_power = np.zeros((time_length, n_days))


    # Creating a figure with multiple subplots, with two rows (one for each type of day)
    fig, ax = plt.subplots(2, 1, sharex = False, sharey = False, figsize = figsize)
    
    # suptitle = 'Results in {}'.format(month)
    # fig.suptitle(suptitle, fontsize = fontsize_title, fontweight = 'bold')
    fig.subplots_adjust(left = 0.1, bottom = 0.1, right = 0.9, top = 0.85, wspace = None, hspace = 0.3)


    for day in days:

        dd = days[day][0]

        pv_production[:, dd]  = pv_production_month[:, mm]  
        ue_consumption[:, dd] = ue_consumption_month_day[:, mm, dd]

        pv_available[:, dd] = pv_production[:, dd] - ue_consumption[:, dd]
        net_load[pv_available[:, dd] < 0, dd] = -pv_available[pv_available[:, dd] < 0, dd]
        pv_available[pv_available[:, dd] < 0, dd]= 0
        
        grid_feed[:, dd], grid_purchase[:, dd], battery_charge[:, dd], battery_discharge[:, dd], battery_energy[:, dd] = \
        battery_optimization(pv_available[:, dd], net_load[:, dd], battery_capacity, battery_specs, grid_power, time_length, time_sim, dt)
    
        shared_power = np.minimum((pv_production + battery_discharge), (ue_consumption + battery_charge))

        number_of_days = months[month]['days_distr'][day]

        pv_energy[mm][dd] = np.sum(pv_production)*dt*number_of_days
        feed_energy[mm][dd] = np.sum(grid_feed)*dt*number_of_days
        purchase_energy[mm][dd] = np.sum(grid_purchase)*dt*number_of_days
        consumption_energy[mm][dd] = np.sum(ue_consumption)*dt*number_of_days
        battery_charge_energy[mm][dd] = np.sum(battery_charge)*dt*number_of_days
        battery_discharge_energy[mm][dd] = np.sum(battery_discharge)*dt*number_of_days
        shared_energy[mm][dd] = np.sum(shared_power)*dt*number_of_days

        tol = 1e-4
        if (np.any(net_load[:, dd] + grid_feed[:, dd] + battery_charge[:, dd] - pv_available[:, dd] - grid_purchase[:, dd] - battery_discharge[:, dd] > tol)): print('ops, node')
        if (np.any(grid_feed[:, dd] > grid_power + tol) or np.any(grid_feed[:, dd] > grid_power + tol)): print('ops, grid power')
        if (np.any(grid_feed[:, dd] > pv_available[:,dd] + tol)): print('ops, grid feed')
        if (np.any(battery_charge[:, dd] > pv_available[:,dd] + tol)): print('ops, battery_charge')


        ax[dd].plot(time_sim + dt/2, pv_available[:, dd], label = 'pv_avaialble')
        ax[dd].plot(time_sim + dt/2, net_load[:, dd], label = 'net_load')
        ax[dd].plot(time_sim + dt/2, grid_feed[:, dd], label = 'grid_feed')
        ax[dd].plot(time_sim + dt/2, grid_purchase[:, dd], label = 'grid_purchase')
        ax[dd].plot(time_sim + dt/2, battery_charge[:, dd], label = 'battery_charge')
        ax[dd].plot(time_sim + dt/2, battery_discharge[:, dd], label = 'battery_discharge')

        ax[dd].bar(time_sim, shared_power[:, dd], color = 'k', width = dt, align = 'edge', label = 'shared power', alpha = 0.3)
        
        title = '{}, {}'.format(month.capitalize(), day)
        ax[dd].set_title(title, fontsize = fontsize_title)

        axtw = ax[dd].twinx()
        axtw.bar(time_sim, battery_energy[:, dd]/battery_capacity*100, color = 'y', width = dt, align = 'edge', label = 'shared power', alpha = 0.3)
        
    ##
    # Making the figure look properly
    for axi in ax.flatten():
        axi.set_xlabel('Time ({})'.format('h'), fontsize = fontsize_labels)
        axi.set_ylabel('Power ({})'.format('kW'), fontsize = fontsize_labels)
        axi.set_xlim([time_sim[0], time_sim[-1]])
        # Set one tick each hour on the x-axis
        axi.set_xticks(list(time_sim[: : int(dt)]))
        axi.tick_params(axis ='both', labelsize = fontsize_ticks)
        axi.tick_params(axis ='x', labelrotation = 0)
        axi.grid()
        axi.legend(loc = 'upper left', fontsize = fontsize_legend, ncol = 2)
    

plt.show()







# # # # Calcolo potenza autoconsumata istantaneamente (condivisa dalla comunità)
# # # P_lgc=np.zeros( (Nint+1,Nmonth) )

# # print(battery_energy[:,0])


# for i in range(time_length):

#     if ((net_load[i,month] + grid_feed[i,month] + battery_charge[i,month] \
#                         - pv_available[i,month] - grid_purchase[i,month] - battery_discharge[i,month])*dt) > 1e-6: print('oh oh')


# # print(battery_discharge[:, month])
# # print(battery_charge[:, month])
# # print(grid_feed[:, month])
# # print(grid_purchase[:, month])

# fig = plt.figure()
# ax = fig.add_subplot()
# # ax.plot(time_sim, battery_discharge[:, month], 'b--', label = 'P_b_d')
# # ax.plot(time_sim, battery_charge[:, month], 'b', label = 'P_b_c')
# ax.plot(time_sim, grid_purchase[:, month], 'r--', label = 'P_g_p')

# ax.plot(time_sim, ue_consumption[:, month], 'b', label = 'ue')
# ax.plot(time_sim, grid_feed[:, month], 'r', label = 'P_g_f')
# # plt.plot(time_sim[:, month], 'g', label = 'pv_ava')
# # plt.plot(net_load[:, month], 'g--', label = 'net_load')
# plt.legend()

# ax_twin = ax.twinx()
# ax_twin.plot(time_sim, battery_energy[:, month], 'y', label = 'E_b')
# plt.legend()

# plt.show()



# # energy_shared = np.minimum((pv + battery_discharge), (ue_consumption + battery_charge))


# # print(np.sum(grid_feed[:, month]))
# # print(np.sum(energy_shared[:, month]))
# # print(np.sum(energy_shared[:, month])/np.sum(pv[:, month] + battery_discharge[:, month])*100)