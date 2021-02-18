# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 16:09:48 2020

@author: giamm
"""
import numpy as np
import csv
from pathlib import Path
from tabulate import tabulate

import datareader
from levenshtein_distance import Leven_dist_comparison

##############################################################################


# This file is used to let the user enter the input parameters from keyboard.
# When a value is not given by the user, a default value is assigned.


##############################################################################

# The base path is saved in the variable basepath, it is used to move among
# directories to find the files that need to be read.
basepath = Path(__file__).parent

# A /Parameters folder is created in order to store the parameters as .csv files
dirname = 'Parameters'

try: Path.mkdir(basepath / dirname)
except Exception: pass 


## Parameters

# Simulation parameters that can be changed
# n_hh              #number of households (-)
# n_people_avg      #average number of members for each household (-)
# ftg_avg           #average footage of each household (m2)
# location          #geographical location: 'north' | 'centre' | 'south'
# power_max         #maximum power available from the grid (contractual power) (W)
# en_class          #energetic class of the appiances: 'A+++' | 'A++' | 'A+' | 'A' | 'B' | 'C' | 'D'
# toll              #tolerance on the displacement of the appliance's daily time-on, i.e. duration (%)
# devsta            #standard deviation of the appliance's daily time-on, i.e. duration (min)
# dt_aggr           #aggregated data timestep (min) 15 | 30 | 60
# q_max             #quantile for the maximum instantaneous load profile (%)
# q_med             #quantile for the medium instantaneous load profile (%)
# q_min             #quantile for the minimum instantaneous load profile (%)
# time_scale        #time-scale for plotting: 'min' | 'h'
# power_scale       #power-scale for plotting: 'W' | 'kW' | 'MW'
# energy_scale      #energy-scale for plotting: 'kWh' | 'MWh'

## Simulation setup

# Parameters for the simulation setup that can be changed (both for PV and battery)
# sim_type          #type of simulation: 'fixed' size | 'parametric'
# size              #fixed size of the system (active if sim_type == 'fixed') (kW)
# size_min          #minimum size of the system (active if sim_type == 'parametric') (kW)
# size_max          #maximum size of the system (active if sim_type == 'parametric') (kW)

################################################################################################################################################################################

## Creating a method to change the parameters entering their values from keyboard
# All the parameters that can be changed are declared as keys of param_dict, while
# the value for each key contains the type (int, float, str) of the parameter, its
# defualt value and its boundaries/possible values and its unit of measure (the latter
# is just to be written in the .csv file where the parameters will be saved).
# The current values of the parameters are read from a .csv file that has been 
# created previously and given as values in a dictionary. If the file does not exist yet,
# default values are applied.
# The user is then asked if there is any change desired in the parameters' values.

def parameters_input():

    '''
    The method creates a dictionary where the parameters (updated to user's input) for the simulation are stored.
    The updated parameters will be saved in a .csv file.
    '''

    # Creating a dictionary that contains all the parameters, their type, default values, etc.
    param_dict = {
        'n_hh': {'type': int, 'default_val': 2, 'min_val': 1, 'max_val': 10000, 'uom': '(units)'},
        # 'toll': {'type': int, 'default_val': 15., 'min_val': 0., 'max_val': 100, 'uom': '(min)'},
        # 'devsta': {'type': int, 'default_val': 2, 'min_val': 1, 'max_val': 100, 'uom': '(min)'},
        # 'q_max': {'type': int, 'default_val': 85, 'min_val': 1, 'max_val': 100, 'uom': '(%)'},
        # 'q_med': {'type': int, 'default_val': 50, 'min_val': 1, 'max_val': 100, 'uom': '(%)'},
        # 'q_min': {'type': int, 'default_val': 15, 'min_val': 1, 'max_val': 100, 'uom': '(%)'},
        # 'n_people_avg': {'type': float, 'default_val': 2.7, 'min_val': 1., 'max_val': 10., 'uom': '(people/household)'},
        'location': {'type': str, 'default_val': 'north', 'possible_values': ['north', 'south', 'centre'], 'uom': '(/)'},
        'power_max': {'type': float, 'default_val': 3., 'min_val': 1., 'max_val': 10., 'uom': '(kW)'},
        'en_class': {'type': str, 'default_val': 'A+', 'possible_values': ['A+++', 'A++', 'A+', 'A', 'B', 'C', 'D'], 'uom': '(/)'},
        'ftg_avg': {'type': float, 'default_val': 100., 'min_val': 10., 'max_val': 1000., 'uom': '(m2)'},
        'dt_aggr':{'type': int, 'default_val': 60, 'possible_values': [15, 30, 60], 'uom': '(min)'},
        # 'time_scale': {'type': str, 'default_val': 'h', 'possible_values': ['min', 'h'], 'uom': '(/)'},
        # 'power_scale': {'type': str, 'default_val': 'kW', 'possible_values': ['W', 'kW', 'MW'], 'uom': '(/)'},
        # 'energy_scale': {'type': str, 'default_val': 'MWh', 'possible_values': ['kWh', 'MWh'], 'uom': '(/)'},
        }

    # Provididng parameters' description
    param_dict['n_hh']['description'] = 'number of households'
    # param_dict['toll']['description'] = 'tollerance on appliances\' duration'
    # param_dict['devsta']['description'] = 'standard deviation on appliances\' duration'
    # param_dict['q_max']['description'] = 'quantile for the maximum instantaneous load profile'
    # param_dict['q_med']['description'] = 'quantile for the medium instantaneous load profile'
    # param_dict['q_min']['description'] = 'quantile for the minimum instantaneous load profile'
    param_dict['dt_aggr']['description'] = 'time-step for the aggregation'
    # param_dict['n_people_avg']['description'] = 'average number of people per household'
    param_dict['ftg_avg']['description'] = 'average footage of the households'
    param_dict['power_max']['description'] = 'maximum (contractual) power'
    param_dict['location']['description'] = 'location (north - centre - south)'
    param_dict['en_class']['description'] = 'energy class of the appliances (A+++ - D)'
    # param_dict['time_scale']['description'] = 'time-scale for plotting'
    # param_dict['power_scale']['description'] = 'power-scale for plotting'
    # param_dict['energy_scale']['description'] = 'energy-scale for plotting'

    # Creating a list that contains all the parameters names (usefull for the inputs from keyboard)
    param_list = list(param_dict.keys())

    # Creating a list of possible commands that will stop the execution of the code
    stop_commands = ['', 'stop', 'done', 'no', 'none']

    # The current values for the parameters a read from the file parameters.csv. If it does not exist yet
    # default values are assigned to the parameters
    params = datareader.read_param('parameters', ';', dirname)

    if not bool(params):
        for param in param_dict: params[param] = param_dict[param]['default_val']

    # Printing the current values for the parameters
    message = '\nThe parameters for the simulation are currently set as follows\n'
    print(message)

    tab = []
    for param in params:
        row = [param, params[param], param_dict[param]['uom'].strip('() '), param_dict[param]['description']]
        tab.append(row)
    
    print(tabulate(tab, headers=['Parameter', 'Value', 'Unit of measure', 'Description']))

    # Starting the input of new values
    message = '\nWould you like to change any parameter?\nPress \'enter\' to avoid or\nEnter \'ok\' to start changing: '
    start_flag = input(message).strip("',.=\"_ ").lower().replace(' ', '_').replace('-', '_')

    if start_flag in stop_commands: message = '\nNo parameter will be changed.\n'
    else: message = '\nUpper/lower cases, underscores, quotation marks can be disregarded.\nPress \'enter\' to stop at any time.'
    print(message)

    # Starting the procedure for updating the values of the parameters
    while start_flag not in stop_commands:

        # Asking for a command-line input in order to change a parameter's value
        param_change = input('\nTo change a parameter write the whole expression (ex. n hh = 100): ') \
            .strip("\"',. ").lower().replace(' ', '_').replace('-', '_')

        # Exiting the loop if a "stop-command" is given
        if param_change in stop_commands: break

        # Finding the equality sign in the expression entered by the user, in order to
        # divide the parameter's name from the value
        index = param_change.find('=')
        param_name = param_change[:index].strip("=\"',. _").lower().replace(' ', '_').replace('-', '_')
        param_val = param_change[index + 1:].strip("=\"',. _ ").lower().replace(' ', '_').replace('-', '_')

        # Assessing if the parameter's name entered by the user is in the parameter's list;
        # otherwise the Leven_dist_comparison method is used to suggest the closest match 
        if param_name not in param_list:
            options = Leven_dist_comparison(param_name, param_list)
            
            if len(options) == 1: message = 'Maybe you mean {} (rewrite the name to confirm): '.format(options[0])
            else : message = 'No match found, please rewrite the parameter\'s name: '
            
            param_name = input(message).strip("=\"',. ").lower().replace(' ', '_').replace('-', '_')

            if param_name in stop_commands: break

        if param_name in stop_commands: break
        elif param_name not in param_list: print('Sorry no quick fix, try again please.'); continue


        # After identifying the parameter that is going to be changed, the value entered by the user
        # is checked to be consistent with the possible values the parameter can assume
        if param_dict[param_name]['type'] == int:

            while True: 
                if param_val in stop_commands: param_val = param_dict[param_name]['default_val']; break

                try: param_val = int(param_val)
                except: param_val = input('Please, enter an integer value for {}: '.format(param_name)) \
                    .strip("=\"',. ").lower().replace(' ', '_').replace('-', '_'); continue
                
                if 'possible_values' not in param_dict[param_name]:
                    
                    low_lim = param_dict[param_name]['min_val']
                    up_lim = param_dict[param_name]['max_val']

                    if param_val >= low_lim and param_val <= up_lim: break
                    else: param_val = input('Please, enter an integer between {} and {}: '.format(low_lim, up_lim)) \
                        .strip("=\"',. ").lower().replace(' ', '_').replace('-', '_'); continue

        
                else:
                    
                    possible_values = param_dict[param_name]['possible_values']
                
                    if param_val in possible_values: break
                    else: param_val = input('Please, enter an integer in {}: '.format(possible_values)) \
                        .strip("=\"',. ").lower().replace(' ', '_').replace('-', '_'); continue
       
        elif param_dict[param_name]['type'] == float:

            low_lim = param_dict[param_name]['min_val']
            up_lim = param_dict[param_name]['max_val']
        
            while True: 
                if param_val in stop_commands: param_val = param_dict[param_name]['default_val']; break

                try: param_val = float(param_val)
                except: param_val = input('Please, enter a number: ') \
                    .strip("=\"',. ").lower().replace(' ', '_').replace('-', '_'); continue
                
                if param_val >= low_lim and param_val <= up_lim: break
                else: param_val = input('Please, enter a number between {} and {}: '.format(low_lim, up_lim)) \
                    .strip("=\"',. ").lower().replace(' ', '_').replace('-', '_'); continue
   
        elif param_dict[param_name]['type'] == str:

            possible_values = param_dict[param_name]['possible_values']
            
            possible_values_low = []
            for value in possible_values: possible_values_low.append(value.lower())
            
            while True: 
                
                if param_val in stop_commands: param_val = param_dict[param_name]['default_val']; break
        
                if param_val in possible_values_low: param_val = possible_values[possible_values_low.index(param_val)]; break
                else: param_val = input('Please, choose between {}: '.format(possible_values)) \
                    .strip("=\"',. ").lower().replace(' ', '_').replace('-', '_'); continue
            
        # Updating the parameter's value
        params[param_name] = param_val
        print('Done: {} changed to {} {}.'.format(param_name, param_val, param_dict[param_name]['uom'].strip('()')))

    # Storing the parameters (updated) in a .csv file
    filename = 'parameters.csv'
    fpath = basepath / dirname 
    with open(fpath / filename , mode='w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=';', quotechar="'", quoting=csv.QUOTE_NONNUMERIC)

        csv_writer.writerow(['Name', 'Value', 'Unit of measure'])

        for param in params:
            csv_writer.writerow([param , params[param], param_dict[param]['uom']])

    # Returning a dictionary with the updated values for the parameters
    return(params)

################################################################################################################################################################################

## Creating a method to change the parameters entering their values from keyboard
# All the parameters that can be changed are declared as keys of param_dict, while
# the value for each key contains the type (int, float, str) of the parameter, its
# defualt value and its boundaries/possible values and its unit of measure (the latter
# is just to be written in the .csv file where the parameters will be saved).
# The current values of the parameters are read from a .csv file that has been 
# created previously and given as values in a dictionary. If the file does not exist yet,
# default values are applied.
# The user is then asked if there is any change desired in the parameters' values.

def simulation_setup(tech):

    '''
    The method creates a dictionary where the simulation setup (updated to user's input) is stored.
    The updated parameters will be saved in a .csv file.

    Input:
        tech - str, the technology about which the simulation setup can changed
    
    Output:
        save_params - dict, containing the updated parameters for the simulation setup
        size_range - list, containing the range of sizes to be explored
    '''


    ## Parameters that can be changed

    # Creating a dictionary that contains all the parameters, their type, default values, etc.
    param_dict = {
        'sim_type': {'type': str, 'default_val': 'fixed', 'possible_values': ['fixed', 'parametric'], 'uom': '(/)'},
        'size': {'type': float, 'default_val': 2, 'min_val': 0.5, 'max_val': 10000, 'uom': '(kW)'},
        'size_min': {'type': float, 'default_val': 2, 'min_val': 0.5, 'max_val': 10000, 'uom': '(kW)'},
        'size_max': {'type': float, 'default_val': 2., 'min_val': 0.5, 'max_val': 10000, 'uom': '(kW)'},
        'n_sizes': {'type': int, 'default_val': 1, 'min_val': 1, 'max_val': 5, 'uom': '(/)'}
        }

    if tech.strip("',.=\"_ ").lower().replace(' ', '_').replace('-', '_') == 'battery':
        param_dict = {
            'sim_type': {'type': str, 'default_val': 'parametric', 'possible_values': ['fixed', 'parametric'], 'uom': '(/)'},
            'size': {'type': float, 'default_val': 2, 'min_val': 0.5, 'max_val': 10000, 'uom': '(kWh)'},
            'size_min': {'type': float, 'default_val': 1, 'min_val': 0.5, 'max_val': 10000, 'uom': '(kWh)'},
            'size_max': {'type': float, 'default_val': 5., 'min_val': 0.5, 'max_val': 10000, 'uom': '(kWh)'},
            'n_sizes': {'type': int, 'default_val': 5, 'min_val': 1, 'max_val': 5, 'uom': '(/)'}
            }

    # Adding a description to each parameter
    param_dict['sim_type']['description'] = 'Type of simulation for {}: \'fixed\' size or \'parametric\''.format(tech)
    param_dict['size']['description'] = 'Fixed size for {}'.format(tech)
    param_dict['size_min']['description'] = 'Minimum size of the {}'.format(tech)
    param_dict['size_max']['description'] = 'Maximum size of the {}'.format(tech)
    param_dict['n_sizes']['description'] = 'Number of sizes of the {} to be evaluated'.format(tech)

    # The current values for the parameters a read from the file parameters.csv. If it does not exist yet
    # default values are assigned to the parameters
    params = datareader.read_param('{}_simulation_setup'.format(tech.strip("',.=\"_ ").lower().replace(' ', '_')), ';', dirname)

    if not bool(params):
        for param in param_dict: params[param] = param_dict[param]['default_val']

    # Printing the current values for the parameters
    message = '\nThe simulation for the {} is currently set as follows\n'.format(tech)
    print(message)

    # The parameters that are printed depend on the type of simulation
    tab = []
    if params['sim_type'] == 'fixed':
        param = 'sim_type'
        tab.append([param, params[param], param_dict[param]['uom'].strip('() '), param_dict[param]['description']])
        param = 'size'
        tab.append([param, params[param], param_dict[param]['uom'].strip('() '), param_dict[param]['description']])
    
    elif params['sim_type'] == 'parametric':
        param = 'sim_type'
        tab.append([param, params[param], param_dict[param]['uom'].strip('() '), param_dict[param]['description']])
        param = 'size_min'
        tab.append([param, params[param], param_dict[param]['uom'].strip('() '), param_dict[param]['description']])
        param = 'size_max'
        tab.append([param, params[param], param_dict[param]['uom'].strip('() '), param_dict[param]['description']])
        param = 'n_sizes'
        tab.append([param, params[param], param_dict[param]['uom'].strip('() '), param_dict[param]['description']])

    print(tabulate(tab, headers=['Parameter', 'Value', 'Unit of measure', 'Description']))


    ## Paramter's update from keyboard

    # Creating a list of possible commands that will stop the execution of the code
    stop_commands = ['', 'stop', 'done', 'no', 'none']

    # Creating a list of the paramters contained in params
    param_list = list(params.keys())

    # Starting the input of new values
    message = '\nWould you like to change any parameter?\nPress \'enter\' to avoid or\nEnter \'ok\' to start changing: '
    start_flag = input(message).strip("',.=\"_ ").lower().replace(' ', '_').replace('-', '_')

    if start_flag in stop_commands: message = '\nNo parameter will be changed.\n'
    else: message = '\nUpper/lower cases, underscores, quotation marks can be disregarded.\nPress \'enter\' to stop at any time.'
    print(message)

    # Starting the procedure for updating the values of the parameters
    while start_flag not in stop_commands:


        # Asking for a command-line input in order to change a parameter's value
        param_change = input('\nWrite the whole expression (ex. sim type = fixed): ').strip("',.=\"_ ").lower().replace(' ', '_').replace('-', '_')

        # Exiting the loop if a "stop-command" is given
        if param_change in stop_commands: break

        # Finding the equality sign in the expression entered by the user, in order to
        # divide the parameter's name from the value
        index = param_change.find('=')
        param_name = param_change[:index].strip("',.=\"_ ").lower().replace(' ', '_').replace('-', '_')
        param_val = param_change[index + 1:].strip("',.=\"_ ").lower().replace(' ', '_').replace('-', '_')

        # Assessing if the parameter's name entered by the user is in the parameter's list;
        # otherwise the Leven_dist_comparison method is used to suggest the closest match 
        if param_name not in param_list:
            options = Leven_dist_comparison(param_name, param_list)
            
            if len(options) == 1: message = 'Maybe you mean {} (rewrite the name to confirm): '.format(options[0])
            else : message = 'No match found, please rewrite the parameter\'s name: '
            
            param_name = input(message).strip("',.=\"_ ").lower().replace(' ', '_').replace('-', '_')

            if param_name in stop_commands: break

        # If a break command is given, the outer loop still has not been broken
        if param_name in stop_commands: break

        # If the parameter name is still not in the list, a new iteration is started
        elif param_name not in param_list: print('Sorry no quick fix, try again please.'); continue


        ## Type of simulation
        # If the the parameter that is going to be changed is the simulation type, also the fixed size/ size range boundaries
        # must be updated, therefore a proper implementation is needed

        if param_name == 'sim_type':

            possible_values = param_dict[param_name]['possible_values']
            possible_values = [value.strip("=\"',. ").lower().replace(' ', '_').replace('-', '_') for value in possible_values]
            default_value = param_dict[param_name]['default_val']

            # if param_val in stop_commands: param_val = default_value
            # If the updated value is not in the possible values, the user is aske dto re-write it, after giving a
            # suggestion. If the value is still not in the possible values, the dafult value is applied
            if param_val not in possible_values:
                options = Leven_dist_comparison(param_val, possible_values)
                
                if len(options) == 1: message = 'Maybe you mean \'{}\' (rewrite the name to confirm): '.format(options[0])
                else : message = 'Please, rewrite the type of simulation you want to perform: '

                param_val = input(message).strip("',.=\"_ ").lower().replace(' ', '_').replace('-', '_')
        
            if param_val in stop_commands: break
            elif param_val not in possible_values: print('Sorry no quick fix, default value will be assigned.'); param_val = default_value

            # Updating the value for the type of simulation
            params[param_name] = param_val
            print('Done: {} changed to {}.'.format(param_name, param_val))

            # As already mentioned, if one wants to change the simulation type, also the size/range boundaries must be changed
            message = '\nIf you want to change the simulation type, you also need to change the size/size-range boundaries'
            print(message)


            ## Sizes range boundaries

            sim_type = param_val

            # For fixed size type, minimum and maximum size are the same
            if sim_type == 'fixed':
                message = 'Enter the size of the {} {}: '.format(tech, param_dict['size']['uom'])
                size = input(message)

                try: default_value = params['size']
                except: default_value = param_dict['size']['default_val']
                low_lim = param_dict['size']['min_val']
                up_lim = param_dict['size']['max_val']
                
                while True: 

                    size = size.strip("',.=\"_ ").lower().replace(' ', '_').replace('-', '_')

                    if size in stop_commands: size = default_value; break

                    try: size = float(size)
                    except: size = input('Please, enter a number: '); continue
                    
                    if size >= low_lim and size <= up_lim: break
                    else: size  = input('Please, enter a number between {} and {}: '.format(low_lim, up_lim)); continue

                params['size'] = size

            # For parametric type, both minimum and maximum size are to be specified
            if sim_type == 'parametric':

                # Minimum size
                message = 'Enter the minimum size of the {} {}: '.format(tech, param_dict['size_min']['uom'])
                size_min = input(message)

                try: default_value = params['size_min']
                except: default_value = param_dict['size_min']['default_val']
                low_lim = param_dict['size_min']['min_val']
                up_lim = param_dict['size_min']['max_val']
                
                while True: 

                    size_min = size_min.strip("',.=\"_ ").lower().replace(' ', '_').replace('-', '_')

                    if size_min in stop_commands: size_min = default_value; break

                    try: size_min = float(size_min)
                    except: size_min = input('Please, enter a number: '); continue
                    
                    if size_min >= low_lim and size_min <= up_lim: break
                    else: size_min  = input('Please, enter a number between {} and {}: '.format(low_lim, up_lim)); continue

                # Maximum size
                message = 'Enter the maximum size of the {} {}: '.format(tech, param_dict['size_max']['uom'])
                size_max = input(message)

                try: default_value = params['size_max']
                except: default_value = param_dict['size_max']['default_val']
                low_lim = param_dict['size_max']['min_val']
                up_lim = param_dict['size_max']['max_val']
                
                while True: 

                    size_max = size_max.strip("',.=\"_ ").lower().replace(' ', '_').replace('-', '_')

                    if size_max in stop_commands: size_max = default_value; break

                    try: size_max = float(size_max)
                    except: size_max = input('Please, enter a number: '); continue
                    
                    if size_max >= low_lim and size_max <= up_lim: break
                    else: size_max  = input('Please, enter a number between {} and {}: '.format(low_lim, up_lim)); continue

                # Making sure size_min is smaller than size_max
                if size_min > size_max: size_min, size_max = size_max, size_min

                # Number of sizes
                message = 'Enter the number of sizes of the {} to be evaluated: '.format(tech)
                n_sizes = input(message).strip("',.=\"_ ").lower().replace(' ', '_').replace('-', '_')

                try: default_value = params['n_sizes']
                except: default_value = param_dict['n_sizes']['default_val']
                low_lim = param_dict['n_sizes']['min_val']
                up_lim = param_dict['n_sizes']['max_val']
                
                while True: 

                    if n_sizes in stop_commands: n_sizes = default_value; break

                    try: n_sizes = int(n_sizes)
                    except: n_sizes = input('Please, enter an integer: ') \
                        .strip("',.=\"_ ").lower().replace(' ', '_').replace('-', '_'); continue
                    
                    if n_sizes >= low_lim and n_sizes <= up_lim: break
                    else: n_sizes  = input('I am good but I can be slow! Please enter an integer between {} and {}: '.format(low_lim, up_lim)) \
                        .strip("',.=\"_ ").lower().replace(' ', '_').replace('-', '_'); continue

                # Storing the updated parameters
                params['size_min'] = size_min
                params['size_max'] = size_max
                params['n_sizes'] = n_sizes
            
            # If the simulation type has been changed and the size/range boundaries have been changed as well there is
            # nothing left to change, so the loop can be broken
            break

        # If the parameter to be changed was not the simulation type, the usual procedure can be followed 
        # to check is the specified value is consisent with the type of the paramter and the minimum/maximum values
        else:

            low_lim = param_dict[param_name]['min_val']
            up_lim = param_dict[param_name]['max_val']
        
            while True: 

                param_val = param_val
                if param_val in stop_commands: param_val = param_dict[param_name]['default_val']; break

                try: param_val = float(param_val)
                except: param_val = input('Please, enter a number: ') \
                    .strip("',.=\"_ ").lower().replace(' ', '_').replace('-', '_'); continue
                
                if param_val >= low_lim and param_val <= up_lim: break
                else: param_val = input('Please, enter a number between {} and {}: '.format(low_lim, up_lim)) \
                    .strip("',.=\"_ ").lower().replace('-', '_'); continue

            params[param_name] = param_val
            print('Done: {} changed to {} {}.'.format(param_name, param_val, param_dict[param_name]['uom'].strip('()')))


    ## Storing the updated parameters
    # Only the ones of interest are stored

    if params['sim_type'] == 'parametric':

        if params['size_min'] == params['size_max'] or params['n_sizes'] == 1:
            save_params = {
                'sim_type': 'fixed',
                'size': params['size_min']
            }
        
        else:
            save_params = {
                'sim_type': params['sim_type'],
                'size_min': params['size_min'],
                'size_max': params['size_max'],
                'n_sizes': int(params['n_sizes']),
            }

    elif params['sim_type'] == 'fixed':
        save_params = {
            'sim_type': params['sim_type'],
            'size': params['size'],
            }


    # Storing the parameters (updated) in a .csv file
    filename = '{}_simulation_setup.csv'.format(tech.strip("',.=\"_ ").lower().replace(' ', '_'))
    fpath = basepath / dirname 
    with open(fpath / filename , mode='w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=';', quotechar="'", quoting=csv.QUOTE_NONNUMERIC)

        csv_writer.writerow(['Name', 'Value'])

        for param in save_params:
            csv_writer.writerow([param , save_params[param]])


    ## Creating the sizes range
    # Once that the simulation setup has been defined, the range for the size can be built and returned to main

    # Updating the simulation type to the user's input
    sim_type = save_params['sim_type']

    # Updating the range's boundaries to the user's input
    if sim_type == 'fixed': size_min, size_max, n_sizes = save_params['size'], save_params['size'], 1
    if sim_type == 'parametric': size_min, size_max, n_sizes = save_params['size_min'], save_params['size_max'], save_params['n_sizes']

    # The boundaries are rounded up to .5 precision
    size_min = int(size_min*2)/2
    size_max = int(size_max*2)/2
    n_sizes = int(n_sizes)

    # The length of the range is evaluated
    size_range_length = max(size_max - size_min, 0.5)

    # Making sure that the difference in size is at least of 0.5 (in case this does not happen, the number of sizes is decreased)
    n_sizes = min(n_sizes, int(2*size_range_length + 1))

    # Creating the size range and making sure that all sizes are rounded up to .5
    size_range = np.linspace(size_min, size_max, n_sizes)
    size_range = [int(size*2)/2 for size in size_range]


#     # The length of the range is evaluated
#     size_range_length = size_max - size_min

#     # The step for the size is chosen depending on the length of the range
#     if size_range_length <= 2.5: d_size = 0.5
#     elif size_range_length > 2.5 and size_range_length <=5: d_size = 1
#     elif size_range_length > 5 and size_range_length <= 10: d_size = 2
#     else: d_size = int(size_range_length/5)

#     # The range is created 
#     size_range = np.arange(size_min, size_max + d_size, d_size)
#     if size_range[-1] != size_max: size_range[-1] = size_max

    return(save_params, list(size_range))
