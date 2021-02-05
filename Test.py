import datareader
import numpy as np
import math
import random
from scipy.interpolate import interp1d
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt
from load_profiler import load_profiler
from load_profile_aggregator_trapz import aggregator
from pathlib import Path
import csv
from tictoc import tic, toc
# from aggregate_load_profiler import aggregate_load_profiler as agr


# The base path is saved in the variable basepath, it is used to move among
# directories to find the files that need to be read.
basepath = Path(__file__).parent

# # Specifying which quantities (load profiles) are going to be stored and plotted
# load_profiles_types = {
#     0: 'Total',
#     1: 'Average',
#     2: 'Maximum',
#     3: 'Medium',
#     4: 'Minimum',
#     }

# header = ['{} load [W]'.format(load_profiles_types[lp]) for lp in load_profiles_types]
# print(['Time (min)'] + header)

#####################################################################################################################ù#
 
# apps, apps_ID, apps_attr = datareader.read_appliances('eltdome_report.csv', ';', 'Input')
# location_dict = {
#         'north': apps_attr['distribution_north'] - (len(apps_attr) - np.size(apps, 1)), 
#         'centre': apps_attr['distribution_centre'] - (len(apps_attr) - np.size(apps, 1)), 
#         'south': apps_attr['distribution_south'] - (len(apps_attr) - np.size(apps, 1)),
#         }

# location = 'south'

# n_hh = 17
# distribution_factors = apps[:, location_dict[location]]

# energies = np.array(17*[100])

# print((n_hh * distribution_factors))

# filename = 'loadprof_ac_wde_sawpOLD.csv'
# data = datareader.read_general(filename, ';', 'Input')
# time = data[:, 0]
# power = data[:, 1]

# print(np.trapz(power, time)*92/1000)

#####################################################################################################################################

# apps, apps_ID, apps_attr = datareader.read_appliances('eltdome_report.csv', ';', 'Input')
# energies = np.random.randint(100, size = (17,))

# apps_classes = {}

# ii = 0
# for app in apps_ID:
#         apps_class = apps_ID[app][apps_attr['class']]
#         if  apps_class not in apps_classes: 
#                 apps_classes[apps_class] = {'id_number': ii}
#                 ii += 1

#         if 'app_list' not in apps_classes[apps_class]: apps_classes[apps_class]['app_list'] = []

#         apps_classes[apps_class]['app_list'].append(apps_ID[app][apps_attr['id_number']])

# n_seasons = 4
# n_hh = 100

# energy_stor = np.random.randint(100, size = (n_seasons, len(apps_ID), n_hh))
# energies_tot_season = np.transpose(np.sum(energy_stor, axis = 2))
# energies_class = np.zeros((len(apps_classes), n_seasons))
# for apps_class in apps_classes:
#         apps_list = apps_classes[apps_class]['app_list']
#         energies_class[apps_classes[apps_class]['id_number'], :] = np.sum(energies_tot_season[apps_list, :], axis = 0)

# print(energies_class)

# print(apps_ID)


#########################################################################################################################################

# myArr1 = np.array([0,1,2,3,4])
# myArr2 = np.array([2,3,5,8,9,45,678,45,34,21,45])

# myMatx = np.array([[0,1,2,3,4], 
#                   [56,78,90,87,4]])

# myMatx[0, -1] = 100 #myArr1
# myMatx[1, -1] = 205 #myArr2

# print(myMatx)

################################################################################################################################################

# apps, apps_ID, apps_attr = datareader.read_appliances('eltdome_report.csv', ';', 'Input')

# time = 1440
# dt = 1
# time_sim = np.arange(0, time, dt)

# # app = 'lighting'

# # # app_wbe (weekly behavior), different usage of the appliance in each type of days: 'wde'|'we','wd'
# # app_wbe = apps_ID[app][apps_attr['week_behaviour']] 

# # # app_sbe (seasonal behavior), different usage of the appliance in each season: 'sawp'|'s','w','ap'
# # app_sbe = apps_ID[app][apps_attr['season_behaviour']] 

# # print('app_wbe: {}, type {}'.format(app_wbe, type(app_wbe)))
# # print('app_sbe: {}, type {}'.format(app_sbe, type(app_sbe)))


# apps_avg_lps = {}
# apps_dcs = {}
# for app in apps_ID:

#         # app_nickname is a 2 or 3 characters string identifying the appliance
#         app_nickname = apps_ID[app][apps_attr['nickname']] 

#         # app_type depends from the work cycle for the appliance: 'continuous'|'no_duty_cycle'|'duty_cycle'|
#         app_type = apps_ID[app][apps_attr['type']]

#         # app_wbe (weekly behavior), different usage of the appliance in each type of days: 'wde'|'we','wd'
#         app_wbe = apps_ID[app][apps_attr['week_behaviour']] 

#         # app_sbe (seasonal behavior), different usage of the appliance in each season: 'sawp'|'s','w','ap'
#         app_sbe = apps_ID[app][apps_attr['season_behaviour']] 

#         # Building the name of the file to be opened and read
#         fname_nickname = app_nickname
#         fname_type = 'avg_loadprof'

#         apps_avg_lps[app] = {}

#         for season in app_sbe:
#                 fname_season = season

#                 for day in app_wbe:
#                         fname_day = day

#                         filename = '{}_{}_{}_{}.csv'.format(fname_type, fname_nickname, fname_day, fname_season)
                        
#                         # Reading the time and power vectors for the load profile
#                         data_lp = datareader.read_general(filename,';','Input')

#                         # Time is stored in hours and converted to minutes
#                         time_lp = data_lp[:, 0] 
#                         time_lp = time_lp*60 

#                         # Power is already stored in Watts, it corresponds to the load profile
#                         power_lp = data_lp[:, 1] 
#                         load_profile = power_lp

#                         # Interpolating the load profile if it has a different time-resolution
#                         if (time_lp[-1] - time_lp[0])/(np.size(time_lp) - 1) != dt: 
#                                 load_profile = np.interp(time_sim, time_lp, power_lp, period = time)

#                         apps_avg_lps[app][(season, day)] = load_profile

        
#         if app_type == 'duty_cycle':
#                 fname_type = 'dutycycle'
#                 filename = '{}_{}.csv'.format(fname_type, fname_nickname)
                
#                 # Reading the time and power vectors for the duty cycle 
#                 data_dc = datareader.read_general(filename, ';', 'Input')

#                 # Time is already stored in  minutes
#                 time_dc = data_dc[:, 0] 

#                 # Power is already stored in Watts, it corresponds to the duty cycle
#                 power_dc = data_dc[:, 1] 
#                 duty_cycle = power_dc
                
#                 # Interpolating the duty-cycle, if it has a different time resolution
#                 if (time_dc[-1] - time_dc[0])/(np.size(time_dc) - 1) != dt:
#                         time_dc = np.arange(time_dc[0], time_dc[-1] + dt, dt)
#                         duty_cycle = np.interp(time_dc, power_dc)

#                 apps_dcs[app] = {'time_dc': time_dc,
#                                  'power_dc': duty_cycle}



######################################################################################################################################################      



        

# seasons = {
#     'winter': (0, 'w'),
#     'spring': (1, 'ap'),
#     'summer': (2, 's'),
#     'autumn': (3, 'ap')
#     }

# days = {
#     'week-day': (0, 'wd'),
#     'weekend-day': (1, 'we')
#     }

# months = {
#     'january': {'id': (0, 'jan'), 'season': 'winter', 'days_distr': {'week-day': 23, 'weekend-day': 8}},
#     'february': {'id': (1, 'feb'), 'season': 'winter', 'days_distr': {'week-day': 20, 'weekend-day': 8}},
#     'march': {'id': (2, 'mar'), 'season': 'winter', 'days_distr': {'week-day': 22, 'weekend-day': 9}},
#     'april': {'id': (3, 'apr'), 'season': 'spring', 'days_distr': {'week-day': 21, 'weekend-day': 9}},
#     'may': {'id': (4, 'may'), 'season': 'spring', 'days_distr': {'week-day': 23, 'weekend-day': 8}},
#     'june': {'id': (5, 'jun'), 'season': 'spring', 'days_distr': {'week-day': 21, 'weekend-day': 9}},
#     'july': {'id': (6, 'jul'), 'season': 'summer', 'days_distr': {'week-day': 22, 'weekend-day': 9}},
#     'august': {'id': (7, 'aug'), 'season': 'summer', 'days_distr': {'week-day': 23, 'weekend-day': 8}},
#     'september': {'id': (8, 'sep'), 'season': 'summer', 'days_distr': {'week-day': 20, 'weekend-day': 10}},
#     'october': {'id': (9, 'oct'), 'season': 'autumn', 'days_distr': {'week-day': 23, 'weekend-day': 8}},
#     'november': {'id': (10, 'nov'), 'season': 'autumn', 'days_distr': {'week-day': 22, 'weekend-day': 8}},
#     'december': {'id': (11, 'dec'), 'season': 'autumn', 'days_distr': {'week-day': 21, 'weekend-day': 10}},
#     }

# seasons_id_list = [seasons[season][0] for season in seasons]
# months_id_list = [months[month]['id'][0] for month in months]

# print(seasons_id_list)
# print(months_id_list)

# days_distr = {}
# for month in months:
#         season = months[month]['season']
#         print(season)
#         if season not in days_distr: days_distr[season] = {'week-day': months[month]['days_distr']['week-day'],
#                                                            'weekend-day': months[month]['days_distr']['weekend-day']}

#         else: 
#                 days_distr[season]['week-day'] += months[month]['days_distr']['week-day']
#                 days_distr[season]['weekend-day'] += months[month]['days_distr']['weekend-day']

# print(days_distr)

      
######################################################################################################################################


# time_length = 24
# n_days = 2
# n_months = 12

# pv_production_month = np.random.randint(100, size = (time_length, n_months))
# ue_consumption_month_day = np.random.randint(100, size = (time_length, n_months, n_days))


# pv_production = np.zeros((time_length, n_days)) 
# ue_consumption = np.zeros((time_length, n_days)) 
# net_load = np.zeros((time_length, n_days)) 
# pv_available = np.zeros((time_length, n_days))
# battery_charge= np.zeros((time_length, n_days))
# battery_discharge = np.zeros((time_length, n_days))
# grid_feed = np.zeros((time_length, n_days))
# grid_purchase = np.zeros((time_length, n_days))
# battery_energy = np.zeros((time_length, n_days))

# mm = 0
# for dd in range(n_days):

#         pv_production[:, dd]  = pv_production_month[:, mm]  
#         ue_consumption[:, dd] = ue_consumption_month_day[:, mm, dd]

#         pv_available[:, dd] = pv_production[:, dd] - ue_consumption[:, dd]

#         print(pv_available[:, dd])


#         net_load[pv_available[:, dd] < 0, dd] = -pv_available[pv_available[:, dd] < 0, dd]
#         pv_available[pv_available[:, dd] < 0, dd]= 0

#         print(pv_available[:, dd])
#         print(net_load[:, dd])
        
#         print('\n\n\n\n')


############################################################################################################################################

# x_min = int(float(input('Min: '))*2)/2

# x_max = int(float(input('Max: '))*2)/2


# x_range_length = x_max - x_min

# if x_range_length <= 2.5: dx = 0.5
# elif x_range_length > 2.5 and x_range_length <=5: dx = 1
# elif x_range_length > 5 and x_range_length <= 10: dx = 2
# else: dx = int(x_range_length/5)

# x_range = np.arange(x_min, x_max + dx, dx)

# print(x_range)
# if x_range[-1] != x_max: x_range[-1] = x_max

# print(x_range)
# print('{} : {} : {}'.format(x_min, dx, x_max))
  

############################################################################################################################################


# possible_values = ['fixed SIZE', 'parametric']
# print(possible_values)

# possible_values = [value.strip("\"',. ").lower().replace(' ', '_') for value in possible_values]
# print(possible_values)

# x = 3; print('x: {}'.format(x))
# y = x; print('y: {}'.format(y))

# x = 4; print('x: {}'.format(x))

# string = '     ___  gggggggghhhhhhh      ..'
# print(string)
# print('I want "{}" all over me'.format(string.strip("\"',. _")))
# print('I want "{}" all over me'.format(string.replace(' ', '_').strip("\"',. ").lower()))

# string = 'n people-avg'
# print(string)
# print(string.replace(' ', '_').replace('-', '_'))

# batt_specs = datareader.read_param('battery_specs.csv', ';', 'Input')
# params = datareader.read_param('parameters.csv', ';', 'Parameters')

# print(batt_specs)
# print(params)

# data_pv = datareader.read_general('pv_production_unit.csv', ';', 'Input')

# time_pv = data_pv[:, 0]
# pv_production_unit = data_pv[:, 1:]

# print(np.shape(time_pv))
# print(np.shape(pv_production_unit))

# time = 24
# dt = 0.50

# time_sim = np.arange(0, time, dt)
# print(time_sim)

# f = interp1d(time_pv, pv_production_unit, kind = 'linear', axis = 0, fill_value = 'extrapolate')

# pv_production_new = f(time_sim)

# print(np.shape(pv_production_new))

# plt.plot(time_pv, pv_production_unit[:, 0], 'b')
# plt.plot(time_sim, pv_production_new[:, 0], 'r--')
# plt.show()



# consumption_seasons = agr_hlp(params)/1000

# months = np.arange(12,1)
# seasons = np.arange(12,3)

# f = interp1d(seasons, consumption_seasons, kind = 'linear', axis = 0, fill_value = 'extrapolate')

# shasha = np.zeros((2,))
# shasha[0] = 'time'
# shasha[1] = 'bubu{}'.format(shasha[0])

# print(shasha)

# mylist = []
# mylist.append('sasha')
# mylist += 'bernie'

# print(mylist)

# aaa = np.array([[1,2,3,4],[1,4,5,6]])
# print(aaa)

# print(aaa[:, 0,np.newaxis])
# aaa = np.concatenate((aaa, aaa[:, 0,np.newaxis]), axis = 1)
# print(aaa)

# dt = 1 #original timestep (min)
# time = 12 #total simulation time (min)
# time_sim = np.arange(0, time, dt)
# time_sim = np.append(time_sim, time)

# load_profile = np.arange(0, time, dt)
# load_profile = np.append(load_profile, load_profile[0])

# new_dt = 2
# new_time = np.arange(0, time, new_dt)

# new_load_profile = np.trapz(load_profile[::2], time_sim[::2])/(new_dt/dt)
# print(new_load_profile)

# plt.plot(time_sim, load_profile)
# plt.bar(new_time, new_load_profile, width = new_dt, fill = False, edgecolor = 'k')
# plt.show()

# print(time_sim)



# ncols = 3
# for i in range(15):
#     print('{}: ({},{})'.format(i, int(i/ncols), i%ncols))
    
# print(int(0.9))
# print(int(0.9999999))

# n_sizes_main = 4
# n_sizes_lead = 3


# data1 = np.random.randint(60, size = (n_sizes_main, n_sizes_lead))
# data2 = np.random.randint(100, size = (n_sizes_main, n_sizes_lead))
# # data1 = np.random.randint(100, size = (n_sizes_main, n_sizes_lead))

# data = np.stack((data1, data2), axis = 2)
# print(np.shape(data))

# fig, ax = plt.subplots(2,2)

# ax[0,0].plot([0,1],[1,1])
# ax[0,1].plot([0,1],[10,12])
# plt.subplot(2,1,2)
# plt.plot([0,2],[12,25])

# plt.show()

# # aaa = np.random.randint(100, size = (2,3,4))
# aaa = np.array([[[1,2,3,4], [4,5,6,7], [7,8,9,0]],
#                 [[4,5,30,43], [344,5,6,7], [6787,8,9,0]],
#                  ])

# print(aaa[:,:,0])

# bbb = np.transpose(aaa, axes = (1, 0, 2))

# # print(np.shape(aaa))
# # print(np.shape(bbb))

# print(bbb[:,:,0])

# plt.plot([0, 0], [1, 5], color = 'b', marker = '', linestyle = '-')
# plt.show()


# ymin_left, ymax_left = 0, 0
# ymin_right, ymax_right = 0, 0


# colors = [(230, 25, 75),
#         (60, 180, 75),
#         (255, 225, 25),
#         (0, 130, 200),
#         (245, 130, 48),
#         (145, 30, 180),
#         (70, 240, 240),
#         (240, 50, 230),
#         (210, 245, 60),
#         (250, 190, 212),
#         (0, 128, 128),
#         (220, 190, 255),
#         (170, 110, 40),
#         (255, 250, 200),
#         (128, 0, 0),
#         (170, 255, 195),
#         (128, 128, 0),
#         (255, 215, 180),
#         (0, 0, 128),
#         (128, 128, 128)]

# # Transforming into rgb triplets
# colors_rgb = []
# for color in colors:
#     color_rgb = []
#     for value in color:
#         color_rgb.append(float(value)/255) 
#     colors_rgb.append(tuple(color_rgb))


# for i_color in range(len(colors_rgb)):
#     color = colors_rgb[i_color]
#     # plt.plot([0, 1], [i_color + 1, i_color + 2], label = '{}'.format(i_color))
#     plt.bar(i_color, 10, align = 'edge', label = '{}'.format(i_color))
#     plt.xlim([0, len(colors_rgb) + 1])

# plt.legend(loc = 'right')

# plt.show()


# from tictoc import tic, toc

# tic()
# k = 0
# for i in range(int(1e6)):
#     k += i
#     k += i + 1

# print('{}: {} s'.format(k, toc()))

# tic()
# k = 0
# for i in range(int(1e6)):
#     k += i
    
# for i in range(int(1e6)):
#     k += i + 1

# print('{}: {} s'.format(k, toc()))

# aa = list(range(0, 10, 1))
# bb = list(range(0, 2, 1))

# print(aa)
# print(bb)

# aa, bb = bb, aa

# print(aa)
# print(bb)



# Creating an /Output folder, if not already existing
dirname = 'Output'

# Creating an /Output/Files folder, if not already existing
subdirname = 'Files'

# Creating a subfolder, if not already existing
subsubdirname = '{}_{}_{}'.format('north', 'D', 10)


# try:

#     data_wd = datareader.read_general('consumption_profiles_month_wd.csv', ';', '/'.join((dirname, subdirname, subsubdirname)))
#     data_we = datareader.read_general('consumption_profiles_month_we.csv', ';', '/'.join((dirname, subdirname, subsubdirname)))


#     consumption_month_wd = data_wd[:, 1:]
#     consumption_month_we = data_we[:, 1:]

#     consumption_month_day = np.stack((consumption_month_wd, consumption_month_we), axis = 2)
#     print(np.shape(consumption_month_day))

#     message = '\nSome load profiles have already been evaluated for the current configuration, do you want to use them?\
#            \nPress \'enter\' to skip and re-evaluate the load profiles \
#            \nEnter \'ok\' to use the available ones: '
#     load_profiler_flag = input(message).strip('\'",._- ').lower()

#     if load_profiler_flag == '': load_profiler_flag = 0
#     else: load_profiler_flag = 1

# except:
#     message = '\nEvaluation of the load profiles for the aggregate of households.'
#     print(message)

# aa = [1]

# bb = [2,3,6,78,900,4]

# print(aa+ bb)

# size_min = 10
# size_max = 10
# n_sizes = 1

# size_min = int(size_min*2)/2
# size_max = int(size_max*2)/2
# n_sizes = int(n_sizes)

# size_range_length = max(size_max - size_min, 0.5)

# # if size_range_length/(n_sizes - 1) < 0.5: n_sizes = int(2*size_range_length + 1)
# n_sizes = min(n_sizes, int(2*size_range_length + 1))
# # print(size_range_length/(n_sizes - 1))

# size_range = np.linspace(size_min, size_max, n_sizes)
# size_range = [int(size*2)/2 for size in size_range]
# print(size_range)

# energy_month = np.random.randint(100, size = (6, 2))

# energy_month = np.array(energy_month, dtype = float)

# energy_month[[3,4],0] = np.nan
# energy_month[1,1] = np.nan

# print(energy_month)
# print(np.shape(energy_month))

# if np.any(np.isnan(energy_month)):
#     print('uuh')
#     index = np.where(np.isnan(energy_month))
#     print(index[0])
#     print(energy_month[index[0]])
#     print([np.isnan(energy_month)])


# y = energy_month.copy()



# nans, x = np.isnan(y), lambda z: z.nonzero()[0]


# print(energy_month)
# print(y)



# print(nans)
# print(x)

# print(y[nans])
# print(x(nans))
# print(x(~nans))
# print(y[~nans])

# y[nans]= np.interp(x(nans), x(~nans), y[~nans])

# plt.plot(y[:,0], 'b')
# plt.plot(energy_month[:,0], 'r--')

# plt.plot(y[:,1], 'c')
# plt.plot(energy_month[:,1], 'm--')

# plt.show()


# y = energy_month.flatten()

# print(y)

# nans, x = np.isnan(y), lambda z: z.nonzero()[0]

# y[nans]= np.interp(x(nans), x(~nans), y[~nans])

# print(y)


# energy_month = y.reshape(6,2)

# print(energy_month)

# energy_month = np.random.randint(100, size = (12, 2))
# energy_month = np.array(energy_month, dtype = float)

# energy_month[[1,2], 0] = np.nan
# energy_month[[4,2], 1] = np.nan
# print(energy_month)

# months = {
#     'january': {'id': (0, 'jan'), 'season': 'winter', 'days_distr': {'week-day': 23, 'weekend-day': 8}},
#     'february': {'id': (1, 'feb'), 'season': 'winter', 'days_distr': {'week-day': 20, 'weekend-day': 8}},
#     'march': {'id': (2, 'mar'), 'season': 'winter', 'days_distr': {'week-day': 22, 'weekend-day': 9}},
#     'april': {'id': (3, 'apr'), 'season': 'spring', 'days_distr': {'week-day': 21, 'weekend-day': 9}},
#     'may': {'id': (4, 'may'), 'season': 'spring', 'days_distr': {'week-day': 23, 'weekend-day': 8}},
#     'june': {'id': (5, 'jun'), 'season': 'spring', 'days_distr': {'week-day': 21, 'weekend-day': 9}},
#     'july': {'id': (6, 'jul'), 'season': 'summer', 'days_distr': {'week-day': 22, 'weekend-day': 9}},
#     'august': {'id': (7, 'aug'), 'season': 'summer', 'days_distr': {'week-day': 23, 'weekend-day': 8}},
#     'september': {'id': (8, 'sep'), 'season': 'summer', 'days_distr': {'week-day': 20, 'weekend-day': 10}},
#     'october': {'id': (9, 'oct'), 'season': 'autumn', 'days_distr': {'week-day': 23, 'weekend-day': 8}},
#     'november': {'id': (10, 'nov'), 'season': 'autumn', 'days_distr': {'week-day': 22, 'weekend-day': 8}},
#     'december': {'id': (11, 'dec'), 'season': 'autumn', 'days_distr': {'week-day': 21, 'weekend-day': 10}},
#     }

# days = {
#     'week-day': (0, 'wd'),
#     'weekend-day': (1, 'we')
#     }

# energy_month = np.zeros((len(months)))

# for month in months:
#     imonth = months[month]['id'][0]

#     weekday_nan_flag = 0

#     for day in days:
#         iday = days[day][0]
#         if iday == 0: weekday = day; #yesterday_days = months[month]['days_distr'][day]
#         elif iday == 1: weekend = day; #tomorrow_days = months[month]['days_distr'][day]
    
#     # print(yesterday)
#     # print(tomorrow)

#     for day in days:
#         iday = days[day][0]

#         energy_this_day = np.random.randint(100)
#         random_prob = np.random.randint(100)
#         if random_prob > 80: energy_this_day = np.nan

#         n_days_today = months[month]['days_distr'][day]

#         if np.isnan(energy_this_day): #print('({},{}) is {}'.format(imonth, iday, energy_month[imonth][iday]))
#             print('{}, {} energy is nan!'.format(month, day))
            
#             if day == weekend:
#                 n_days_weekday = months[month]['days_distr'][weekday]
#                 energy_this_day = energy_month[imonth]/n_days_weekday*n_days_today

#             if day == weekday:
#                 weekday_nan_flag = 1
        
#         elif day == weekend and weekday_nan_flag == 1 :

#             # yesterday = days.keys()[days.values().index(0)] 
#             n_days_weekday = months[month]['days_distr'][weekday]
#             energy_month[imonth] = energy_this_day/n_days_today*n_days_weekday

#         energy_month[imonth] += energy_this_day

# print(energy_month.T)

# if np.any(np.isnan(energy_month)):

#     nans, fnan = np.isnan(energy_month), lambda z: z.nonzero()[0]
#     energy_month[nans]= np.interp(fnan(nans), fnan(~nans), energy_month[~nans], period = np.size(energy_month))

# print(energy_month.T)


# energy_month = np.zeros((len(months)))

# for month in months:
#     imonth = months[month]['id'][0]

#     weekday_nan_flag = 0

#     for day in days:
#         iday = days[day][0]
#         if iday == 0: weekday = day; #yesterday_days = months[month]['days_distr'][day]
#         elif iday == 1: weekend = day; #tomorrow_days = months[month]['days_distr'][day]
    
#     # print(yesterday)
#     # print(tomorrow)

#     for day in days:
#         iday = days[day][0]

#         power_this_day = np.random.randint(100, size = (24,))
        
#         random_prob = np.random.randint(100)
#         if random_prob > 80: power_this_day = np.nan

#         n_days_today = months[month]['days_distr'][day]
#         energy_month[imonth] += np.nansum(power_this_day)*n_days_today

#         # print('{}, {} energy is {}'.format(month, day, np.nansum(power_this_day)*n_days_today))
    

#         if np.any(np.isnan(power_this_day)): 

#             print('{}, {} energy is nan'.format(month, day))
            
#             if day == weekday:
#                 weekday_nan_flag = 1

#             elif day == weekend and weekday_nan_flag == 0:
#                 n_days_weekday = months[month]['days_distr'][weekday]
#                 energy_month[imonth] += energy_month[imonth]/n_days_weekday*n_days_today

#             elif day == weekend and weekday_nan_flag == 1:
#                 energy_month[imonth] = np.nan

        
#         elif day == weekend and weekday_nan_flag == 1 :

#             # yesterday = days.keys()[days.values().index(0)] 
#             n_days_weekday = months[month]['days_distr'][weekday]
#             energy_month[imonth] += energy_month[imonth]/n_days_today*n_days_weekday

#         # print(energy_month[imonth])
#         # energy_month[imonth] += energy_this_day

# print(energy_month.T)

# if np.any(np.isnan(energy_month)):

#     nans, fnan = np.isnan(energy_month), lambda z: z.nonzero()[0]
#     energy_month[nans]= np.interp(fnan(nans), fnan(~nans), energy_month[~nans], period = np.size(energy_month))

# print(energy_month.T)


# aaa = np.zeros((15,))
# bbb = np.zeros((15,))
# ccc = np.zeros((15,))
# ddd = np.ones((15,))

# aaa[:] = np.nan
# # bbb[:] = np.nan
# # ccc[:] = np.nan
# # ddd[:] = np.nan

# eee = np.minimum(aaa + bbb, ccc + ddd)
# print(eee)

###################################################################################################################################################### HIGH NUMBER LOAD PROFILES

# # Time-step, total time and vector of time from 00:00 to 23:59 (one day) (min)
# dt = 1 
# time = 1440 
# time_sim = np.arange(0,time,dt) 

# # Creating a dictionary to be passed to the various methods, containing the time discretization
# time_dict = {
#     'time': time,
#     'dt': dt,
#     'time_sim': time_sim,
#     }

# # Uploading apps attributes 
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

        
#     if app_type == 'duty_cycle':
#         fname_type = 'dutycycle'
#         filename = '{}_{}.csv'.format(fname_type, fname_nickname)
        
#         # Reading the time and power vectors for the duty cycle 
#         data_dc = datareader.read_general(filename, ';', 'Input')

#         # Time is already stored in  minutes
#         time_dc = data_dc[:, 0] 

#         # Power is already stored in Watts, it corresponds to the duty cycle
#         power_dc = data_dc[:, 1] 
#         duty_cycle = power_dc
        
#         # Interpolating the duty-cycle, if it has a different time resolution
#         if (time_dc[-1] - time_dc[0])/(np.size(time_dc) - 1) != dt:
#                 time_dc = np.arange(time_dc[0], time_dc[-1] + dt, dt)
#                 duty_cycle = np.interp(time_dc, power_dc)

#         apps_dcs[app] = {'time_dc': time_dc,
#                         'duty_cycle': duty_cycle}


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
#     'toll': 5,
#     'devsta': 0,
#     'ftg_avg': 100,
# }

# # # app = 'vacuum_cleaner'
# # # app = 'air_conditioner'
# # app = 'electric_oven'
# # # app = 'microwave_oven'
# # # app = 'fridge'
# # app = 'freezer'
# # app = 'washing_machine'
# app = 'dish_washer'
# # # app = 'tumble_drier'
# # app = 'electric_boiler'
# # # app = 'hifi_stereo'
# # app = 'dvd_reader'
# # app = 'tv'
# # # app = 'iron'
# # # app = 'pc'
# # # app = 'laptop'
# # # app = 'lighting'

# day = 'we'
# season = 's'

# # app_wbe (weekly behavior), different usage of the appliance in each type of days: 'wde'|'we','wd'
# app_wbe = apps_ID[app][apps_attr['week_behaviour']] 

# # app_sbe (seasonal behavior), different usage of the appliance in each season: 'sawp'|'s','w','ap'
# app_sbe = apps_ID[app][apps_attr['season_behaviour']] 

# key_season = 'sawp' 
# if len(app_sbe) > 1: key_season = season

# # Default choice (no different behaviour for different types of day):
# # if the appliance has got different profiles in different days of the week, this will be changed
# key_day = 'wde' 
# if len(app_wbe) > 1: key_day = day

# avg_load_profile = apps_avg_lps[app][(key_season, key_day)]

# app_ID = apps_ID[app][apps_attr['id_number']]
# T_on = apps[app_ID, apps_attr['time_on'] - (len(apps_attr) - np.size(apps, 1))]
# indices = int(T_on/2/dt)

# load_profile = np.zeros((np.size(time_sim)))

# for i in range(int(5e4)):
#     load_profile += load_profiler(time_dict, app, day, season, appliances_data, **params)

# plt.plot(time_sim, load_profile/np.max(load_profile), 'b', linewidth = 2)
# plt.bar(time_sim, avg_load_profile/np.max(avg_load_profile), width = dt, color = 'm', alpha = 0.5)
# plt.plot(time_sim, np.roll(load_profile/np.max(load_profile), +indices), 'r', linewidth = 1)
# # plt.ylim(0.65, 1.01)
# plt.show()



################################################################################################################################################### Interpolation/aggregtion

# dt = 1
# time = 24
# time_sim = np.arange(0, time, dt)
# # power = np.random.randint(100, size = (100,))



# mm = 7

# data_pv = datareader.read_general('pv_production_unit.csv', ';', 'Input')
# time_pv = data_pv[:, 0]
# pv_production_unit = data_pv[:, 1:]

# pv_size = 12

# pv_production = pv_production_unit[:, mm]*pv_size

# location = 'north'
# en_class = 'A+'
# n_hh = 12

# # Checking if there is already a file where the load profiles for this configuration have been stored
# dirname = 'Output'
# subdirname = 'Files'
# subsubdirname = '{}_{}_{}'.format(location, en_class, n_hh)

# data_wd = datareader.read_general('consumption_profiles_month_wd.csv', ';', '/'.join((dirname, subdirname, subsubdirname)))
# # data_we = datareader.read_general('consumption_profiles_month_we.csv', ';', '/'.join((dirname, subdirname, subsubdirname)))
# consumption_month_wd = data_wd[:, 1:]
# # consumption_month_we = data_we[:, 1:]

# consumption = consumption_month_wd[:, mm]

# consumption_shared = np.random.rand(np.size(time_sim))*(np.max(consumption)/10)*0

# plt.figure()
# plt.plot(time_sim, pv_production)
# plt.plot(time_sim, consumption)
# plt.plot(time_sim, consumption_shared)
# plt.show()

# net_production = pv_production - consumption_shared
# net_production[net_production < 0] = 0

# power_available = net_production - consumption
# power_available[power_available < 0] = 0

# net_load = consumption + consumption_shared - pv_production


# power_available2 = pv_production - (consumption_shared + consumption)
# power_available2[power_available2 < 0] = 0
# # print(power_available[power_available > 0])
# # print(-net_load[power_available > 0])


# plt.figure()
# # plt.plot(time_sim + dt/2, pv_production, 'bs-', label = 'production')
# # plt.plot(time_sim + dt/2, consumption + consumption_shared, 'gs-', 'consumption')
# plt.plot(time_sim + dt/2, net_production, 'bd--', label = 'net_production')
# plt.plot(time_sim + dt/2, power_available, 'rd--', label = 'power available')
# plt.plot(time_sim + dt/2, power_available2, 'rs-.', label = 'power_available2')
# plt.plot(time_sim + dt/2, net_load, 'gd-.', label = 'net_load')
# plt.legend()
# plt.show()





# dt_aggr = 1
# time_aggr = np.arange(0, time, dt_aggr)
# power_intp = np.interp(time_aggr, time_sim, power, period = time)
# # power_aggr = aggregator({'dt': dt, 'time': time, 'time_sim': time_sim}, power, dt_aggr)




# plt.bar(time_sim, power, color = 'b', align = 'edge', alpha = 0.5, width = dt)
# plt.bar(time_aggr, power_intp, color = 'r', align = 'edge', alpha = 0.5, width = dt_aggr)
# # plt.bar(time_aggr, power_aggr, align = 'edge', color = 'm', alpha = 0.5, width = dt_aggr)

# print(np.trapz(np.append(power, power[0]), np.append(time_sim, time)))
# print(np.sum(power_intp)*dt_aggr)
# # print(np.sum(power_aggr)*dt_aggr)

# plt.show()




##################################################################################################################################################

# location = 'north'
# n_hh = 12
# en_class = 'A+'



# seasons = {
#     'winter': (0, 'w'),
#     'spring': (1, 'ap'),
#     'summer': (2, 's'),
#     'autumn': (3, 'ap')
#     }

# days = {
#     'week-day': (0, 'wd'),
#     'weekend-day': (1, 'we')
#     }

# # While the routine that evaluates the load profiles works with typical profiles for each season, 
# # the optimization procedure works with tyipcal profiles for each month.  
# # A reference year is considered, in which the first day (01/01) is a monday. 
# # Therefore, conventionally considering that winter lasts from january to march 
# # spring from april to june, summer from july to september and autumn from october
# # to december, each month has got the following number of weekdays and weekend days.

# # Creating a dictionary for the months
# months = {
#     'january': {'id': (0, 'jan'), 'season': 'winter', 'days_distr': {'week-day': 23, 'weekend-day': 8}},
#     'february': {'id': (1, 'feb'), 'season': 'winter', 'days_distr': {'week-day': 20, 'weekend-day': 8}},
#     'march': {'id': (2, 'mar'), 'season': 'winter', 'days_distr': {'week-day': 22, 'weekend-day': 9}},
#     'april': {'id': (3, 'apr'), 'season': 'spring', 'days_distr': {'week-day': 21, 'weekend-day': 9}},
#     'may': {'id': (4, 'may'), 'season': 'spring', 'days_distr': {'week-day': 23, 'weekend-day': 8}},
#     'june': {'id': (5, 'jun'), 'season': 'spring', 'days_distr': {'week-day': 21, 'weekend-day': 9}},
#     'july': {'id': (6, 'jul'), 'season': 'summer', 'days_distr': {'week-day': 22, 'weekend-day': 9}},
#     'august': {'id': (7, 'aug'), 'season': 'summer', 'days_distr': {'week-day': 23, 'weekend-day': 8}},
#     'september': {'id': (8, 'sep'), 'season': 'summer', 'days_distr': {'week-day': 20, 'weekend-day': 10}},
#     'october': {'id': (9, 'oct'), 'season': 'autumn', 'days_distr': {'week-day': 23, 'weekend-day': 8}},
#     'november': {'id': (10, 'nov'), 'season': 'autumn', 'days_distr': {'week-day': 22, 'weekend-day': 8}},
#     'december': {'id': (11, 'dec'), 'season': 'autumn', 'days_distr': {'week-day': 21, 'weekend-day': 10}},
#     }

# # The days distribution in the seasons can be evaluated as well
# days_distr = {}
# for month in months:
#         season = months[month]['season']
        
#         if season not in days_distr: days_distr[season] = {'week-day': months[month]['days_distr']['week-day'],
#                                                            'weekend-day': months[month]['days_distr']['weekend-day']}

#         else: 
#                 days_distr[season]['week-day'] += months[month]['days_distr']['week-day']
#                 days_distr[season]['weekend-day'] += months[month]['days_distr']['weekend-day']

# # Storing some useful quantities
# n_months = len(months)
# n_seasons = len(seasons)
# n_days = len(days)

# # Creating two lists to properly slice the consumption data when interpolating from seasonal to monthly values
# months_slice = np.arange(0, n_months)
# seasons_slice = np.arange(0, n_months, int(n_months/n_seasons))

# # Storing all these dictionaries in a new dict that is passed to the various methods
# auxiliary_dict = {
#     'seasons': seasons,
#     'n_seasons': n_seasons,
#     'months': months,
#     'n_months': n_months,
#     'days': days,
#     'n_days': n_days,
#     'days_distr': days_distr,
#     }




# ## Time discretization
# # Time is discretized in steps of one hour (according to the definition of shared energy in the Decree Law 162/2019 - "Mille proroghe")

# # Total time of simulation (h) - for each typical day
# time = 24

# # Timestep for the simulation (h) (the keyboard input dt_aggr is given in minutes)
# dt = dt_aggr/60

# # Vector of time, from 00:00 to 23:59, i.e. 24 h
# time_sim = np.arange(0, time, dt)
# time_length = np.size(time_sim)

# # Storing all the elements needed for the time-discretization in a dictionary that is passed to the various methods
# time_dict = {
#     'time': time,
#     'dt': dt,
#     'time_sim': time_sim,
#     'time_length': time_length,
#     }



# ### Input data

# # Maximum power from the grid (total for the aggregate of households) (kW)
# grid_power_max = power_max*n_hh


# ## Battery specification
# # The battery specifications are stored in a file that can be read using the method read_param
# # from the module datareader.py, that will return a dictionary

# battery_specs = datareader.read_param('battery_specs.csv', ';', 'Input')

# # Storing the information about the various technologies considered in a dictionary
# # that is passed to the various methods
# technologies_dict = {
#     'grid_power_max': grid_power_max,
#     'battery_specs': battery_specs,
#     }


# ## Unit production from the photovoltaic installation (kWh/h/kWp)
# # The unit production during each hour (kWh/h/kWp) from the photovoltaic installation can 
# # be read using the method read_general from the module datareader.py, that returns a 2d-array
# # containing the time vector on the first columns and the unit production during each hour in each
# # month in the other columns

# data_pv = datareader.read_general('pv_production_unit.csv', ';', 'Input')

# time_pv = data_pv[:, 0]
# pv_production_unit = data_pv[:, 1:]







# ## Consumption from the aggregate of households 

# # Checking if there is already a file where the load profiles for this configuration have been stored
# dirname = 'Output'
# subdirname = 'Files'
# subsubdirname = '{}_{}_{}'.format(location, en_class, n_hh)

# # If the files exist, the user is asked if it is needed to re-evaluate the load profiles or the available ones can be used
# try:

#     data_wd = datareader.read_general('consumption_profiles_month_wd.csv', ';', '/'.join((dirname, subdirname, subsubdirname)))
#     data_we = datareader.read_general('consumption_profiles_month_we.csv', ';', '/'.join((dirname, subdirname, subsubdirname)))

#     consumption_month_wd = data_wd[:, 1:]
#     consumption_month_we = data_we[:, 1:]

#     consumption_month_day = np.stack((consumption_month_wd, consumption_month_we), axis = 2)




######################################################################################################################################################### METHOD FOR READING PV DATA
# filename = 'PVGIS_Data'
# dirname = 'Input'
# delimit = ','

# dirname = dirname.strip()

# filename = filename.strip()
# if not filename.endswith('.csv'): filename = filename + '.csv'

# fpath = basepath / dirname 

# # month_list = []
# # hour_list = []
# # power_list = []

# pv_production = np.zeros((24, 12))
# pv_production_count_days = np.zeros((24, 12))

# tic()

# # try:
# if True:
#     with open(fpath / filename, mode = 'r') as csv_file:
#         csv_reader = csv.reader(csv_file, delimiter = delimit)

        
#         row_before = []
#         headers_flag = 1
#         for row in csv_reader:

            
#             if not row == [] and row[0][0].isdigit():

        
#                 if headers_flag == 1:
#                     headers = row_before
#                     headers = [header.strip().lower().replace(' ', '_') for header in headers]
#                     headers_flag = 0
                    

#                 time_row = row[headers.index('time')]
#                 month = int(time_row[4:6]) - 1
#                 hour = int(time_row[9:11])
#                 power = float(row[headers.index('p')])/1000

#                 pv_production[hour, month] += power
#                 pv_production_count_days[hour, month] += 1
#                 # date_list.append([time_row[0:4], time_row[4:6], time_row[6:8] , time_row[9:11]])
#                 # month_list.append(time_row[4:6])
#                 # hour_list.append(time_row[9:11])
#                 # power_list.append(row[headers.index('p')])
              
            
#             if headers_flag == 1: row_before = row 

# # print(headers_flag)      
# # except:
    
# #     print('Unable to open the file')
# #     # print('Unable to open this file')

# # print(date_list[100])
# # print(date_list[200])

# # print(power_list[100])
# # print(power_list[200])
# # Creating a 2D-array containing the data(time in the first column and power in the second one)
# # month = np.array(month_list, dtype = int)
# # hour = np.array(hour_list, dtype = int)
# # power = np.array(power_list, dtype = float)

# # print(type(date[0, 0]))
# # print(type(date[0, -1]))
# # print(type(power[0]))

# # print(date[0, :])
# # print(power[0])



# pv_production = np.column_stack((np.arange(0, 24, 1), pv_production/pv_production_count_days))

# print(toc())


# # Storing the parameters (updated) in a .csv file
# filename = 'pv_production_try2.csv'
# fpath = basepath / dirname 
# with open(fpath / filename , mode='w', newline='') as csv_file:
#     csv_writer = csv.writer(csv_file, delimiter = ';', quotechar="'", quoting = csv.QUOTE_NONNUMERIC)

#     csv_writer.writerow(['Time (h)'] + ['Month {} (W)'.format(i) for i in range(12)])

#     for row in pv_production:
#         csv_writer.writerow(row)


# for i_month in range(np.size(pv_production, axis = 1)):
#     print(i_month)

#     for i_hour in range(np.size(pv_production, axis = 0)):
#         print(i_hour)

#         print(month[month == i_month])

#         pv_production[hour, month] = np.average(power[np.all(month == i_month and hour == i_hour)])

# # print(pv_production)

