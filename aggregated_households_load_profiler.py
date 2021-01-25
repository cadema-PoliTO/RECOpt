# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 10:53:11 2020

@author: giamm
"""

import numpy as np
# import matplotlib.pyplot as plt
import random
# import math
import csv
from pathlib import Path

import datareader
import parameters_input as inp
import plot_generator  as plot
from house_load_profiler import house_load_profiler as hlp
from load_profile_aggregator import aggregator as agr
from tictoc import tic, toc

###############################################################################


# This is the _main_ file of a routine that generates the electric load-profile
# for an aggregate of a certain number of households. The total simulation time 
# is 1440 min, while the resolution (timestep) is 1 min.
# The load profile and the daily energy consumption are evaluated for a number
# of households and appliances, according to the availability of each appliance
# in each household (depending on a distribution factor and on the geographical
# location that is chosen). The load profile from all the households is then 
# aggregated and the result is shown with a different time resolution (dt_aggr).
# This is done for each season, both for weekdays and weekend days. The energy 
# consumption from each appliance during the whole year is also evaluated.


###############################################################################

# The basepath of the file is stored in a variable 
basepath = Path(__file__).parent


## Parameters needed in the simulation

# Simulation parameters that can be changed by keyboard input
# n_hh              #number of households (-)
# n_people_avg      #average number of members for each household (-)
# ftg_avg           #average footage of each household (m2)
# location          #geographical location: 'north' | 'centre' | 'south'
# power_max         #maximum power available from the grid (contractual power) (W)
# en_class          #energetic class of the appiances: 'A+++' | 'A++' | 'A+' | 'A' | 'B' | 'C' | 'D'
# dt_aggr           #aggregated data timestep (min) 15 | 30 | 60

# Parameters that are to be changed manually
# toll              #tolerance on the displacement of the appliance's daily time-on, i.e. duration (%)
# devsta            #standard deviation of the appliance's daily time-on, i.e. duration (min)
# q_max             #quantile for the maximum instantaneous load profile (%)
# q_med             #quantile for the medium instantaneous load profile (%)
# q_min             #quantile for the minimum instantaneous load profile (%)
# time_scale        #time-scale for plotting: 'min' | 'h'
# power_scale       #power-scale for plotting: 'W' | 'kW' | 'MW'
# energy_scale      #energy-scale for plotting: 'kWh' | 'MWh'

def aggregate_load_profiler():


    ## Parameters

    # Updating the parameters according to the keyboard input by calling the parameters_input() method
    params = inp.parameters_input()

    # Some more "insiders" paramters that are to be changed manually
    params['devsta'] = 2
    params['toll'] = 15 
    params['q_max'] = 80
    params['q_med'] = 50
    params['q_min'] = 20

    # Updating the parameters' values (only for those parameters that are needed here)

    # Number of households considered (-)
    n_hh = params['n_hh'] 

    # Geographical location: 'north' | 'centre' | 'south'
    location = params['location']

    # Energetic class of the appiances: 'A+++' | 'A++' | 'A+' | 'A' | 'B' | 'C' | 'D'
    en_class = params['en_class']

    # Time-step used to aggregated the results (min): 1 | 5 | 10 | 15 | 10 | 30 | 45 | 60
    dt_aggr = params['dt_aggr']

    # Quantile for the evaluation of maximum, medium and minimum power demands at each time-step
    q_max = params['q_max']
    q_med = params['q_med']
    q_min = params['q_min']




    #### Starting the simulation
    tic()


    ## Time 
    # Time discretization for the simulation    

    # Time-step, total time and vector of time from 00:00 to 23:59 (one day) (min)
    dt = 1 
    time = 1440 
    time_sim = np.arange(0,time,dt) 

    # Time vector for the aggregation of the results (min)
    time_aggr = np.arange(0,time,dt_aggr) 


    ## Input data for the appliances
    # Appliances' attributes, energy consumptions and user's coefficients 

    apps, apps_ID, apps_attr = datareader.read_appliances('eltdome_report.csv', ';', 'Input')
    # apps is a 2d-array in which, for each appliance (rows) and attribute value is given (columns)
    # apps_ID is a dictionary in which, for each appliance (key), its ID number,type,week and seasonal behavior (value)
    # apps_attr is a dictionary in which the name of each attribute (value) is linked to its columns number in apps (key)

    ec_yearly_energy, ec_levels_dict = datareader.read_enclasses('classenerg_report.csv', ';', 'Input')
    # ec_yearly_energy is a 2d-array in which for each appliance, its yearly energy consumption is given for each energetic class
    # ec_levels_dict is a dictionary that links each energetic level (value) to its columns number in ec_yearly_energy

    coeff_matrix, seasons_dict = datareader.read_enclasses('coeff_matrix.csv',';','Input')
    # coeff_matrix is a 2d-array in which for each appliance, its coefficient k, related to user's behaviour in different seasons, is given
    # seasons_dict is a dictionary that links each season (value) to its columns number in coeff_matrix

    # Average daily load profiles (from MICENE, REMODECE) and typical duty cycle diagrams (from CESI)
    # Rather than loading the proper file each time the load_profiler method is called, all the input data are loaded here
    # and stored in two dicionaries (one for the average daily load profiles and one for the typical duty cycle diagrams) 
    # that are passed each time to the functions

    apps_avg_lps = {}
    apps_dcs = {}
    for app in apps_ID:

        # Storing some useful variables needed to identify the file to be loaded and read

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

        # Initializing the dictionary (value) related to the current appliance - app (key)
        apps_avg_lps[app] = {}

        # Running through different seasons if the appliance's usage changes according to the season 
        # if app_sbe == 'sawp' there will only be one iteration
        for season in app_sbe:
            fname_season = season

            # Running through different day-types if the appliance's usage changes according to the day-type
            # if app_wbe == 'wde' there will only be one iteration
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

                # Storing the load profile in the proper element of the dictionary
                apps_avg_lps[app][(season, day)] = load_profile

        # Loading the duty cycle diagram for the appliance of "duty-cycle" type
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

            # Storing time and power vectors of the duty cycle
            apps_dcs[app] = {'time_dc': time_dc,
                            'duty_cycle': duty_cycle}

    # Creating a dictionary to pass such data to the various methods
    appliances_data = {
        'apps': apps,
        'apps_ID': apps_ID,
        'apps_attr': apps_attr,
        'ec_yearly_energy': ec_yearly_energy,
        'ec_levels_dict': ec_levels_dict,
        'coeff_matrix': coeff_matrix,
        'seasons_dict': seasons_dict,
        'apps_avg_lps': apps_avg_lps,
        'apps_dcs': apps_dcs,
        }


    ## Building the appliances availability matrix
    # A 2d-array is built in which, for each household (columns) it is
    # shown which appliances are available for the household, according to the 
    # distribution factor of each appliance in a given location (1 if available, 0 otherwise)

    # Initializing the array
    apps_availability = np.zeros((len(apps_ID),n_hh))
    number_of_apps = np.zeros((len(apps_ID)))

    # A dictionary that relates each location (key) to the related columns in apps(for distribution factors)
    location_dict = {
        'north': apps_attr['distribution_north'] - (len(apps_attr) - np.size(apps, 1)), 
        'centre': apps_attr['distribution_centre'] - (len(apps_attr) - np.size(apps, 1)), 
        'south': apps_attr['distribution_south'] - (len(apps_attr) - np.size(apps, 1)),
        }

    # Building the matrix
    for app in apps_ID:

        # The ID number of the appliance is stored in a variable since it will be used man times
        app_ID = apps_ID[app][apps_attr['id_number']]  
        
        # Extracting the distribution factor for the appliance in the current geographical location
        distr_fact = apps[app_ID, location_dict[location]] 

        # Evaluating the number of households in which the appliance is available      
        n_apps_app_type = int(np.round(n_hh*distr_fact))
        number_of_apps[app_ID] = n_apps_app_type

        # Extracting randomly the households in which the appliance is available, from the total number of households  
        samp = random.sample(list(range(0, n_hh)), n_apps_app_type) 
        
        # Assigning the appliance's availability to the households present in the sample
        apps_availability[apps_ID[app][0],samp] = 1 


    ## Building seasons and week dictionaries 
    # This is done in order to explore all the seasons and, for each season two types of days (weekday and weekend)      
    seasons = {'winter': (0, 'w'), 'spring': (1, 'ap'), 'summer': (2, 's'), 'autumn': (3, 'ap')}
    days = {'week-day': (0, 'wd'), 'weekend-day': (1, 'we')}

    #  A reference year is considered, in which the first day (01/01) is a monday. 
    # Therefore, conventionally considering that winter lasts from 21/12 to 20/03, 
    # spring from 21/03 to 20/06, summer from 21/06 to 20/09 and winter from 21/09
    # to 20/12, each seasons has got the following number of weekdays and weekend days.
    days_distr = {'winter': {'week-day': 64, 'weekend-day': 26},
                'spring': {'week-day': 66, 'weekend-day': 26},
                'summer': {'week-day': 66, 'weekend-day': 26},
                'autumn': {'week-day': 65, 'weekend-day': 26}
                }



    ### Evaluating the load profiles for all the households

    # The aggregated load profiles are evaluated for both a week day and a weekend
    # day for the four seasons and the seasonal energy consumption from each 
    # appliance, for each household are evaluated.

    # First, some quantities are initialized, that will be useful for storing the results

    # Quantile are evaluated for the load profiles. It means that for each timestep
    # the maximum power (demanded by less than 15% of the households), the median
    # power (demanded by less than 50% of the households) and the minimum (demanded 7
    # by less than 85% of the households).
    nmax = int(np.round(n_hh*q_max/100))
    nmed = int(np.round(n_hh*q_med/100))
    nmin = int(np.round(n_hh*q_min/100)) 

    # A random sample of n_samp houses is also extracted in order to plot, for each 
    # of them, the load profile during one day for each season, for each day-type.
    n_samp = 5
    if n_hh < n_samp: n_samp = n_hh

    # A random sample is extracted from the total number of households
    samp = random.sample(list(range(0, n_hh)), n_samp) 

    # A list where to store the header for a .csv file is initialized
    sample_lp_header = samp

    # Storing some useful quantities into variables
    n_seasons = len(seasons)
    n_days = len(days)
    n_apps = len(apps_ID)
    n_time_sim = np.size(time_sim)
    n_time_aggr = np.size(time_aggr)

    # Specifying which quantities (load profiles) are going to be stored and plotted
    load_profiles_types = {
        0: 'Total',
        1: 'Average',
        2: 'Maximum',
        3: 'Medium',
        4: 'Minimum',
        }
    lps_header = ['{} load (W)'.format(load_profiles_types[lp]) for lp in load_profiles_types]

    # Creating 3d-arrays where to store each type of load profile, for each season (axis = 0),
    # each time-step (axis = 1) and type of day (axis = 2)
    lp_tot_stor = np.zeros((n_seasons, n_time_aggr, n_days))
    lp_avg_stor = np.zeros((n_seasons, n_time_aggr, n_days))
    lp_max_stor = np.zeros((n_seasons, n_time_aggr, n_days))
    lp_med_stor = np.zeros((n_seasons, n_time_aggr, n_days))
    lp_min_stor = np.zeros((n_seasons, n_time_aggr, n_days))

    # Creating a 3d.array where to store the energy consumption, for each season (axis = 0),
    # each appliance (axis = 1) and household (axis = 2)
    energy_stor = np.zeros((n_seasons, n_apps, n_hh))

    # Creating a 4d-array where to store the load profiles for the random sample households (axis = 0)
    # for each season (axis = 1), each time-step (axis = 2) and type of day (axis = 3)
    lp_samp_stor = np.zeros((n_samp, n_seasons, n_time_sim, n_days))


    ## Evaluating the load profiles for the households, using the methods created, for each season and type of day

    # Running through the seasons using a for loop
    for season in seasons:

        # Season ID number and nickname to be used in the following
        ss = seasons[season][0]
        season_nickname = seasons[season][1]

        # Running through the day-types using a for loop
        for day in days:
            
            # Day-type ID number and nickname to be used in the following
            dd = days[day][0]
            day_nickname = days[day][1]
            

            ## Generating the load profile (lp) for all the households, 
            # According to the availability of each appliance in the household
            
            # A 2d-array is initalized, containing the load profile in time (axis = 0) for each household (axis = 1)
            households_lp = np.zeros((n_time_sim, n_hh)) 

            # A 2d-array is initialized, containing the energy consumed in one day by each appliance (axis = 0) for each household (axis = 1)
            apps_energy = np.zeros((n_apps, n_hh)) 
            
            # The house_load_profiler routine is applied to each household
            for household in range(0, n_hh):
                households_lp[:, household], apps_energy[:, household] = \
                    hlp(apps_availability[:, household], day_nickname, season_nickname, appliances_data, **params)


            ## Aggregating the load profiles and storing the results
        
            # Aggregating the load profiles with a different timestep
            aggr_households_lp = agr(households_lp, dt_aggr)

            # Evaluating the sum of all the load profiles, for each time-step
            aggregate_lp = np.sum(aggr_households_lp, axis = 1) 

            # Sorting the values of the power from all the households in an increasing way, for each time-step
            sorted_lp = np.sort(aggr_households_lp, axis = 1) 
            
            # Evaluating the maximum, medium and minimum load profiles (quantile), for each time-step
            quantile_lp_max = sorted_lp[:, nmax]
            quantile_lp_med = sorted_lp[:, nmed]
            quantile_lp_min = sorted_lp[:, nmin]

            # Evaluating and storing the random sample load profiles for the current season and day-type
            lp_samp_stor[:, ss, :, dd] = households_lp[:, samp].transpose()

            # Saving all the load profiles in the proper array for the current season and day-type
            lp_tot_stor[ss, :, dd] = aggregate_lp
            lp_avg_stor[ss, :, dd] = aggregate_lp/n_hh
            lp_max_stor[ss, :, dd] = quantile_lp_max
            lp_med_stor[ss, :, dd] = quantile_lp_med
            lp_min_stor[ss, :, dd] = quantile_lp_min
            

            ## Evaluating the energy consumed by the appliances 
            
            # For the current season, it is given by the energy consumed in the current day-type
            # multiplied by the number of days of this type present in the season 
            n_days_season = days_distr[season][day]
            energy_stor[ss, :, :] += apps_energy*n_days_season



    ### Saving the results as .csv files

    # Creating an /Output folder, if not already existing
    dirname = 'Output'

    try: Path.mkdir(basepath / dirname)
    except Exception: pass

    # Creating an /Output/Files folder, if not already existing
    subdirname = 'Files'

    try: Path.mkdir(basepath / dirname / subdirname)
    except Exception: pass

    message = '\nThe results are ready and are now being saved in {}/{}.\n'.format(dirname, subdirname)
    print(message)  

    # Creating a subfolder, if not already existing
    subsubdirname = '{}_{}_{}'.format(location, en_class, n_hh)

    try: Path.mkdir(basepath / dirname / subdirname / subsubdirname)
    except Exception: pass


    ## Running through the seasons and the day-types

    # Storing the load profiles and energy consumptions as .csv files
    for season in seasons:
        ss = seasons[season][0]

        for day in days:
            dd = days[day][0]
            day_nickname = days[day][1]

            # Saving the total, average, maximum, medium, minimum (i.e. quantile) load profiles in a .csv file
            filename = '{}_{}_{}_{}_{}_lps_aggr.csv'.format(location, en_class, n_hh, season, day_nickname)
            fpath = basepath / dirname / subdirname / subsubdirname
                
            with open(fpath / filename, mode = 'w', newline = '') as csv_file:
                csv_writer = csv.writer(csv_file, delimiter = ';', quotechar = "'", quoting = csv.QUOTE_NONNUMERIC)
                csv_writer.writerow(['Time (min)'] + lps_header)

                for ii in range(np.size(time_aggr)):
                    csv_writer.writerow([time_aggr[ii], lp_tot_stor[ss, ii, dd], lp_avg_stor[ss, ii, dd], lp_max_stor[ss, ii, dd], lp_med_stor[ss, ii, dd], lp_min_stor[ss, ii, dd]])

            # Saving the random sample load profiles in a file, after giving a different time-step
            filename = '{}_{}_{}_{}_{}_lps_sample.csv'.format(location, en_class, n_hh, season, day_nickname)
            fpath = basepath / dirname / subdirname / subsubdirname

            with open(fpath / filename, mode = 'w', newline = '') as csv_file:
                csv_writer = csv.writer(csv_file, delimiter = ';', quotechar = "'", quoting = csv.QUOTE_NONNUMERIC)
                csv_writer.writerow(sample_lp_header)

                for ii in range(np.size(time_sim)):
                    csv_writer.writerow([time_sim[ii]] + list(lp_samp_stor[:, ss, ii, dd]))
            
        # Saving the energy consumed by the appliances in each household, for each season
        filename = '{}_{}_{}_{}_energy.csv'.format(location, en_class, n_hh, season)
        fpath = basepath / dirname / subdirname / subsubdirname
        
        with open(fpath / filename, mode='w', newline='') as csv_file:
            csv_writer = csv.writer(csv_file, delimiter=';', quotechar="'", quoting=csv.QUOTE_NONNUMERIC)
            csv_writer.writerow(['App_name', 'App_nickname'] + list(range(0,n_hh)))
        
            for app in apps_ID:
                csv_writer.writerow([app,apps_ID[app][1]]+list(energy_stor[ss, apps_ID[app][0], :]))   



    ### Post-processing of the results

    dirname = 'Output'

    # Creating an /Output/Figures folder, if not already existing
    subdirname = 'Figures'

    message = '\nThe results are now being plotted, figures are saved in {}/{}.\n'.format(dirname, subdirname)
    print(message)  

    try: Path.mkdir(basepath / dirname / subdirname)
    except Exception: pass

    # Creating a subfolder, if not already existing
    subsubdirname = '{}_{}_{}'.format(location, en_class, n_hh)

    try: Path.mkdir(basepath / dirname / subdirname / subsubdirname)
    except Exception: pass


    ## Running through the seasons and calling the methods in plot_generator to create figures to be saved

    for season in seasons:

        ss = seasons[season][0]

        # Total load profiles       
        plot_specs = {
            0: ['bar', 'Total'],
            }

        fig_specs = {
            'suptitle': '\n'.join(('\nTotal load profile during one day', \
                        'for {} households with {} energetic class in the {} of Italy'\
                        .format(n_hh, en_class, location.capitalize()))),
        }

        powers = lp_tot_stor[np.newaxis, ss, :, :]
                            
        fig = plot.seasonal_load_profiles(time_aggr, powers, season, plot_specs, fig_specs, appliances_data, **params)

        filename = '{}_{}_{}_{}_aggr_loadprof.png'.format(location, en_class, season, n_hh)
        fpath = basepath / dirname / subdirname / subsubdirname
        
        fig.savefig(fpath / filename) 

        # Average load profile and quantiles
        plot_specs = {
        0: ['plot', 'Average'],
        1: ['bar', 'Max'],
        2: ['bar', 'Med'],
        3: ['bar', 'Min']
        }

        fig_specs = {
            'suptitle': '\n'.join(('\nAverage load profile and quantile during one day', \
                        'for {} households with {} energetic class in the {} of Italy'\
                        .format(n_hh, en_class, location.capitalize()))),
        }

        powers = np.stack((lp_avg_stor[ss, :, :],
                        lp_max_stor[ss, :, :],
                        lp_med_stor[ss, :, :],
                        lp_min_stor[ss, :, :]),
                        axis = 0)

        fig = plot.seasonal_load_profiles(time_aggr, powers, season, plot_specs, fig_specs, appliances_data, **params)

        filename = '{}_{}_{}_{}_avg_quant_loadprof.png'.format(location, en_class, season, n_hh)
        fpath = basepath / dirname / subdirname / subsubdirname
        
        fig.savefig(fpath / filename)

        # Random sample load profiles
        plot_specs = {}
        for ii in range(n_samp):
            plot_specs[ii] = ['plot','Household: {}'.format(samp[ii])]

        fig_specs = {
            'suptitle': '\n'.join(('\nRandom sample load profiles during one day', \
                        'for {} households with {} energetic class in the {} of Italy'\
                        .format(n_hh, en_class, location.capitalize()))),
        }

        powers = lp_samp_stor[:, ss, :, :]

        fig = plot.seasonal_load_profiles(time_sim, powers, season, plot_specs, fig_specs, appliances_data, **params)

        filename = '{}_{}_{}_{}_sample_loadprof.png'.format(location, en_class, season, n_hh)
        fpath = basepath / dirname / subdirname / subsubdirname
        
        fig.savefig(fpath / filename)


    # Total energy consumption by season
    energies_season = np.transpose(np.sum(energy_stor, axis = 2))

    fig_specs = {
        'suptitle': '\n'.join(('Total energy consumption from appliances by season',
                    'for {} households with {} energetic class in the {} of Italy'.format(n_hh, en_class, location.capitalize())))
        }

    fig = plot.seasonal_energy(apps_ID, energies_season, fig_specs, appliances_data, **params)

    filename = '{}_{}_{}_season_tot_energy_apps.png'.format(location, en_class, n_hh)
    fpath = basepath / dirname / subdirname / subsubdirname

    fig.savefig(fpath / filename)


    # Yearly total energy consumption
    fig_specs = {
        'suptitle': '\n'.join(('\nTotal energy consumption from the appliances for one year',
                    'for {} households with {} energetic class in the {} of Italy'.format(n_hh, en_class, location.capitalize()))),
        'text': 'on',
    }

    energies_year = np.sum(energies_season, axis = 1)

    fig = plot.yearly_energy(apps_ID, energies_year, fig_specs, appliances_data, **params)

    filename = '{}_{}_{}_year_tot_energy_apps.png'.format(location, en_class, n_hh)
    fpath = basepath / dirname / subdirname / subsubdirname

    fig.savefig(fpath / filename)


    # Yearly average energy consumption
    fig_specs = {
        'suptitle': '\n'.join(('\nAverage energy consumption from the appliances (considering ownership) for one year',
                    'for {} households with {} energetic class in the {} of Italy'.format(n_hh, en_class, location.capitalize()))),
        'text': 'off',
    }

    # Dividing each energy consumption by the number of units of each appliance
    number_of_apps[number_of_apps == 0] = -1
    avg_energies = energies_year / number_of_apps

    fig = plot.yearly_energy(apps_ID, avg_energies, fig_specs, appliances_data, **params, energy_scale = 'kWh')

    filename = '{}_{}_{}_year_avg_energy_apps.png'.format(location, en_class, n_hh)
    fpath = basepath / dirname / subdirname / subsubdirname

    fig.savefig(fpath / filename)


    # Percentage of energy consumption by classes of appliances
    apps_classes = {}
    apps_classes_labels = {}

    ii = 0
    for app in apps_ID:
        apps_class = apps_ID[app][apps_attr['class']]
        if  apps_class not in apps_classes: 
            apps_classes[apps_class] = {'id_number': ii}
            apps_classes_labels[apps_class] = [ii, ii]
            ii += 1

        if 'app_list' not in apps_classes[apps_class]: apps_classes[apps_class]['app_list'] = []

        apps_classes[apps_class]['app_list'].append(apps_ID[app][apps_attr['id_number']])

    energies_class = np.zeros((len(apps_classes), n_seasons))
    for apps_class in apps_classes:
        apps_list = apps_classes[apps_class]['app_list']
        energies_class[apps_classes[apps_class]['id_number'], :] = np.sum(energies_season[apps_list, :], axis = 0)

    fig_specs = {
        'suptitle': '\n'.join(('\nPercentage of total energy consumption from classes of appliances by season',
                    'for {} households with {} energetic class in the {} of Italy'.format(n_hh,en_class,location.capitalize()))),
    }

    energies_perc = energies_class/np.sum(energies_class)*100
    
    fig = plot.seasonal_energy_pie(apps_classes_labels, energies_perc, fig_specs, appliances_data, **params)
    filename = '{}_{}_{}_season_perc_energy_apps.png'.format(location, en_class, n_hh)
    fpath = basepath / dirname / subdirname / subsubdirname

    fig.savefig(fpath / filename)
        
        


    #### End of the simulation   

    message = '\nEnd.\nTotal time: {0:.3f} .\n'.format(toc())
    print(message)

    return(lp_tot_stor)


pp = aggregate_load_profiler()