# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 10:45:12 2020

@author: giamm
"""

from pathlib import Path
import numpy as np
import csv
import math

##############################################################################


# This scripted is used to create methods that properly read files


##############################################################################

# The base path is saved in the variable basepath, it is used to move among
# directories to find the files that need to be read.
basepath = Path(__file__).parent


##############################################################################
# if True:

#     dirname = 'Parameters'
#     filename = 'parameters'
#     delimit = ';'

def read_param(filename,delimit,dirname):
    
    ''' The function reads from a .csv file in which some parameters are saved. 
    The file has the parameter's name in the first column, its value 
    in the second one and its unit of measure (uom) in the third one.
        
    Inputs:
        filname - string containing the name of the file (extension of the file: .dat)
        delimit - string containing the delimiting element
        dirname - name of the folder where to find the file to be opened and read
        
    Outputs:
        params - dict, containing the parameters (keys) and their values, as entered by 
                 by the user and stored in the .csv file        
    '''
    
    dirname = dirname.strip()

    filename = filename.strip()
    if not filename.endswith('.csv'): filename = filename + '.csv'
    
    fpath = basepath / dirname 
    
    params = {}
    
    try:
        with open(fpath / filename, mode='r') as csv_file:
            csv_reader = csv.reader(csv_file,delimiter=delimit,quotechar="'")
            
            header_row = 1
            for row in csv_reader:

                if header_row == 1:

                    for ii in range(len(row)): row[ii] = row[ii].lower().strip().replace(' ', '_')
                    header = {
                        'name': row.index('name'),
                        'value': row.index('value'),
                    }

                    header_row = 0
                    continue

                else:

                    param_name = row[header['name']]
                    param_val = row[header['value']]
                    
                    try: 
                        param_val = int(param_val)
                    except: 
                        try: param_val = float(param_val)
                        except: param_val = param_val
                        
                    params[param_name] = param_val
                          
    except:
        print('Unable to open this file')
    
    # print('Im returning params: {}'.format(params))
    return(params)

 
##############################################################################
    
def read_general(filename,delimit,dirname):
    
    ''' The function reads from a .csv file in which the header is a single row
        
    
    Inputs:
        filname - string containing the name of the file (extension of the file: .dat)
        delimit - string containing the delimiting element
        dirname - name of the folder where to find the file to be opened and read
        
    Outputs:
        data - 2d-array containing the values in the file'
    '''
    
    dirname = dirname.strip()

    filename = filename.strip()
    if not filename.endswith('.csv'): filename = filename + '.csv'
    
    fpath = basepath / dirname 
    
    data_list=[]
    
    try:
        with open(fpath / filename, mode='r') as csv_file:
            csv_reader = csv.reader(csv_file,delimiter=delimit)
            next(csv_reader, None) 
            for row in csv_reader:
                
                data_list.append(row)              
                
    except:
        
        print('Unable to open this file')
    
    # Creating a 2D-array containing the data(time in the first column and power in the second one)
    data = np.array(data_list,dtype='float')
    return(data)


##############################################################################

def read_appliances(filename, delimit, dirname):
    
    ''' The function reads from a .csv file that contains all the appliances
    
    Inputs:
        filname - string containing the name of the file (extension of the file: .dat)
        delimit - string containing the delimiting element
        dirname - name of the folder where to find the file to be opened and read
        
    Outputs:
        app - 2D-array containing, for each appliance, its attributes' values
        app_ID - dictionary containing for ID (key) the related appliance's name'
        app_attributes - dictionary containing for each attribute (columns in app) its description and unit of measure
    '''
    
    dirname = dirname.strip()

    filename = filename.strip()
    if not filename.endswith('.csv'): filename = filename + '.csv'
    
    fpath = basepath / dirname  
    
    # Initializing a list that contains, for each appliance, the numerical values of its attributes 
    apps_list = [] 

    # Initializing a dictionary that contains, for each appliance, its ID number, nickname, type, weekly and seasonal behaviour and class
    apps_ID = {} 
    
    # Reading the CSV file
    with open(fpath / filename, mode='r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=delimit)
        
        # Initializing a flag (header_row) that is used to properly treat the header
        header_row = 1
        for row in csv_reader:
            if header_row == 1:
                header = row

                # Skipping the second row, that contains the units of measure of the attributes
                next(csv_reader)

                header_row = 0
                continue
            
            else:

                # Storing only the numeircal values in app_list
                apps_list.append(row[7:])

                # Storing the non-numerical values (ID, nickname and so on) in app_ID
                apps_ID[row[1].lower().replace(' ', ';')] = (int(row[0]), row[2], row[3], row[4].split(','), row[5].split(','), row[6])
            
            
    # Creating a dictionary that contains the appliances' attributes
    apps_attributes = {}
    
    ii = 0
    for attr in header:

        if attr.lower().replace(' ',';') == 'name': continue
        apps_attributes[attr.lower().replace(' ', '_')] = ii
        ii += 1
    
    # Creating a 2D-array containing appliances and attributes
    apps = np.array(apps_list, dtype = 'float')
    return(apps, apps_ID, apps_attributes)


##############################################################################

def read_enclasses(filename, delimit, dirname):
    
    ''' The function reads from a .csv file that contains, for each appliance, its yearly energy consumption (kWh/year) for every energetic class
    
    Inputs:
        filname - string containing the name of the file (extension of the file: .dat)
        delimit - string containing the delimiting element
        dirname - name of the folder where to find the file to be opened and read
        
    Outputs:
        enclass_en - 2D-array containing, for each appliance (rows), and for each energetic class (columns) the yearly energy consumption (kWh/year)
        enclass_levels - dictionary containing for each energetic class (columns in app) its level 
    '''
    
    apps_ID = read_appliances('eltdome_report', ';', 'Input')[1]
    
    dirname = dirname.strip()
  
    filename = filename.strip()
    if not filename.endswith('.csv'): filename = filename + '.csv'
    
    fpath = basepath / dirname  
    
     # Initializing a dictionary that contains, for each appliance, its nominal yearly energy consumption for all energy classes
    enclass_dict = {}
        
    # Reading the CSV file
    with open(fpath / filename, mode = 'r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter = delimit)
    
        # Initializing a flag (header_row) that is used to properly treat the header
        header_row = 1
        for row in csv_reader:
            if header_row == 1:
                header = row

                # Skipping the second row, that contains the units of measure of the attributes
                next(csv_reader)

                header_row = 0
                continue
            
            else:
                # Storing the yearly energy consumption for each energy class (values), for each appliance (keys)
                enclass_dict[row[0].lower().replace(' ', ';')] = row[1:]
                            
    # Creating a dictionary for energetic classes' levels
    enclass_levels = {}
    
    ii = 0
    for attr in header[1:]:
        enclass_levels[attr] = ii
        ii += 1

    # Number of energetic classes (needed in the next step)
    enclass_n = len(enclass_levels)

    # A list is initialized, where the yearly energy consumptions for each appliance can be stored,
    # after being sorted is the same order as apps_ID (values are initialized to -1, so that exceptions will 
    # occur later on, if an appliance is present in apps_ID but not in the enclass_dict)
    enclass_sorted = [[-1]*enclass_n]*len(apps_ID)
    
    for app in enclass_dict:
        enclass_sorted[apps_ID[app][0]] = enclass_dict[app]
    
   
    # Creating a 2D-array containing appliances and nominal energy consumptions
    enclass_en = np.array(enclass_sorted, dtype='float') 
    return(enclass_en, enclass_levels)


##############################################################################

def read_energy(filename, delimit, dirname):
    
    ''' The function reads from a .csv file that contains, for each appliance, its yearly energy consumption (kWh/year) for every household.
    
    Inputs:
        filname - string containing the name of the file (extension of the file: .dat)
        delimit - string containing the delimiting element
        dirname - name of the folder where to find the file to be opened and read
        
    Outputs:
        energy - 2d-array containing in each cell the value of the seasonal energy consumption from each appliance (rows) for each household (columns)
    '''
   
    dirname = dirname.strip()
  
    filename = filename.strip()
    if not filename.endswith('.csv'): filename = filename + '.csv'
    
    fpath = basepath / dirname 
    
    energy_list = [] #list containing, for each appliance,its nominal energy consumption 
        
    # Reading the CSV file
    with open(fpath / filename, mode = 'r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter = delimit)
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                line_count += 1
                continue
            
            else:
                energy_list.append(row[2:])
                            
            line_count += 1
    
    
    # Creating a 2D-array containing appliances and nominal energy consumptions
    energy = np.array(energy_list,dtype='float')     
    return(energy)

##############################################################################



# apps, apps_ID, apps_attributes = read_appliances('eltdome_report', ';', 'Input')

# print(apps)
# print(apps_ID)
# print(apps_attributes)

# en_class, en_class_levels = read_enclasses('classenerg_report', ';', 'Input')

# print(en_class)
# print(en_class_levels)

# coeff_matrix, seasons_dict = read_enclasses('coeff_matrix', ';', 'Input')

# print(coeff_matrix)
# print(seasons_dict)