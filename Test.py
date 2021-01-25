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

apps, apps_ID, apps_attr = datareader.read_appliances('eltdome_report.csv', ';', 'Input')

time = 1440
dt = 1
time_sim = np.arange(0, time, dt)

# app = 'lighting'

# # app_wbe (weekly behavior), different usage of the appliance in each type of days: 'wde'|'we','wd'
# app_wbe = apps_ID[app][apps_attr['week_behaviour']] 

# # app_sbe (seasonal behavior), different usage of the appliance in each season: 'sawp'|'s','w','ap'
# app_sbe = apps_ID[app][apps_attr['season_behaviour']] 

# print('app_wbe: {}, type {}'.format(app_wbe, type(app_wbe)))
# print('app_sbe: {}, type {}'.format(app_sbe, type(app_sbe)))


apps_avg_lps = {}
apps_dcs = {}
for app in apps_ID:

        # app_nickname is a 2 or 3 characters string identifying the appliance
        app_nickname = apps_ID[app][apps_attr['nickname']] 

        # app_type depends from the work cycle for the appliance: 'continuous'|'no_duty_cycle'|'duty_cycle'|
        app_type = apps_ID[app][apps_attr['type']]

        # app_wbe (weekly behavior), different usage of the appliance in each type of days: 'wde'|'we','wd'
        app_wbe = apps_ID[app][apps_attr['week_behaviour']] 

        # app_sbe (seasonal behavior), different usage of the appliance in each season: 'sawp'|'s','w','ap'
        app_sbe = apps_ID[app][apps_attr['season_behaviour']] 

        # Building the name of the file to be opened and read
        fname_nickname = app_nickname
        fname_type = 'avg_loadprof'

        apps_avg_lps[app] = {}

        for season in app_sbe:
                fname_season = season

                for day in app_wbe:
                        fname_day = day

                        filename = '{}_{}_{}_{}.csv'.format(fname_type, fname_nickname, fname_day, fname_season)
                        
                        # Reading the time and power vectors for the load profile
                        data_lp = datareader.read_general(filename,';','Input')

                        # Time is stored in hours and converted to minutes
                        time_lp = data_lp[:, 0] 
                        time_lp = time_lp*60 

                        # Power is already stored in Watts, it corresponds to the load profile
                        power_lp = data_lp[:, 1] 
                        load_profile = power_lp

                        # Interpolating the load profile if it has a different time-resolution
                        if (time_lp[-1] - time_lp[0])/(np.size(time_lp) - 1) != dt: 
                                load_profile = np.interp(time_sim, time_lp, power_lp, period = time)

                        apps_avg_lps[app][(season, day)] = load_profile

        
        if app_type == 'duty_cycle':
                fname_type = 'dutycycle'
                filename = '{}_{}.csv'.format(fname_type, fname_nickname)
                
                # Reading the time and power vectors for the duty cycle 
                data_dc = datareader.read_general(filename, ';', 'Input')

                # Time is already stored in  minutes
                time_dc = data_dc[:, 0] 

                # Power is already stored in Watts, it corresponds to the duty cycle
                power_dc = data_dc[:, 1] 
                duty_cycle = power_dc
                
                # Interpolating the duty-cycle, if it has a different time resolution
                if (time_dc[-1] - time_dc[0])/(np.size(time_dc) - 1) != dt:
                        time_dc = np.arange(time_dc[0], time_dc[-1] + dt, dt)
                        duty_cycle = np.interp(time_dc, power_dc)

                apps_dcs[app] = {'time_dc': time_dc,
                                 'power_dc': duty_cycle}



        

# print(len(apps_avg_lps))    
# print(len(apps_avg_lps['lighting']))
# print(len(apps_avg_lps['washing_machine']))
# print(np.shape(apps_avg_lps['tv'][('sawp','wde')]))

# print(len(apps_dcs))
# print(len(apps_dcs['washing_machine']))
# print(len(apps_dcs['dish_washer']))
# print(np.shape(apps_dcs['tumble_drier']['power_dc'])) 


        # # Default choice (no different behaviour for different types of day):
        # # if the appliance has got different profiles in different days of the week, this will be changed
        # fname_day = 'wde' 
        # if len(app_wbe) > 1: fname_day = day

        # # Default choice (no different behaviour for different seasons): 
        # # if the appliance has got different profiles in different seasons, this will be changed
        # fname_season = 'sawp' 
        # if len(app_sbe) > 1: fname_season = season

        




      



  

       


