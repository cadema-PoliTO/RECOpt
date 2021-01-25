import datareader
import numpy as np

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


time_length = 24
n_days = 2
n_months = 12

pv_production_month = np.random.randint(100, size = (time_length, n_months))
ue_consumption_month_day = np.random.randint(100, size = (time_length, n_months, n_days))


pv_production = np.zeros((time_length, n_days)) 
ue_consumption = np.zeros((time_length, n_days)) 
net_load = np.zeros((time_length, n_days)) 
pv_available = np.zeros((time_length, n_days))
battery_charge= np.zeros((time_length, n_days))
battery_discharge = np.zeros((time_length, n_days))
grid_feed = np.zeros((time_length, n_days))
grid_purchase = np.zeros((time_length, n_days))
battery_energy = np.zeros((time_length, n_days))

mm = 0
for dd in range(n_days):

        pv_production[:, dd]  = pv_production_month[:, mm]  
        ue_consumption[:, dd] = ue_consumption_month_day[:, mm, dd]

        pv_available[:, dd] = pv_production[:, dd] - ue_consumption[:, dd]

        print(pv_available[:, dd])


        net_load[pv_available[:, dd] < 0, dd] = -pv_available[pv_available[:, dd] < 0, dd]
        pv_available[pv_available[:, dd] < 0, dd]= 0

        print(pv_available[:, dd])
        print(net_load[:, dd])
        
        print('\n\n\n\n')

  

       


