import datareader
import numpy as np
import math
import random
from scipy.interpolate import interp1d
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt
# from aggregate_load_profiler import aggregate_load_profiler as agr

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

size_min = 10
size_max = 10
n_sizes = 1

size_min = int(size_min*2)/2
size_max = int(size_max*2)/2
n_sizes = int(n_sizes)

size_range_length = max(size_max - size_min, 0.5)

# if size_range_length/(n_sizes - 1) < 0.5: n_sizes = int(2*size_range_length + 1)
n_sizes = min(n_sizes, int(2*size_range_length + 1))
# print(size_range_length/(n_sizes - 1))

size_range = np.linspace(size_min, size_max, n_sizes)
size_range = [int(size*2)/2 for size in size_range]
print(size_range)