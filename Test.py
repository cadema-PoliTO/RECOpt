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

apps, apps_ID, apps_attr = datareader.read_appliances('eltdome_report.csv', ';', 'Input')
energies = np.random.randint(100, size = (17,))

apps_classes = {}

ii = 0
for app in apps_ID:
        apps_class = apps_ID[app][apps_attr['class']]
        if  apps_class not in apps_classes: 
                apps_classes[apps_class] = {'id_number': ii}
                ii += 1

        if 'app_list' not in apps_classes[apps_class]: apps_classes[apps_class]['app_list'] = []

        apps_classes[apps_class]['app_list'].append(apps_ID[app][apps_attr['id_number']])

n_seasons = 4
n_hh = 100

energy_stor = np.random.randint(100, size = (n_seasons, len(apps_ID), n_hh))
energies_tot_season = np.transpose(np.sum(energy_stor, axis = 2))
energies_class = np.zeros((len(apps_classes), n_seasons))
for apps_class in apps_classes:
        apps_list = apps_classes[apps_class]['app_list']
        energies_class[apps_classes[apps_class]['id_number'], :] = np.sum(energies_tot_season[apps_list, :], axis = 0)

print(energies_class)

print(apps_ID)