# -*- coding: utf-8 -*-
"""
Created on Fri Feb 05 10:24:12 2021

@author: giamm
"""
import numpy as np
import csv
from pathlib import Path

import datareader

##############################################################################


# This scripted is used to read specifically a .csv file containing hourly
# PV production data for different years, as downloaded from PVGIS.
# The data are elaborated and hourly production profiles for 12 typical days
# (one for each month) in a year are stored


##############################################################################

# The base path is saved in the variable basepath, it is used to move among
# directories to find the files that need to be read.
basepath = Path(__file__).parent


## Original data file

# Filename
filename = 'PVGIS_Data'
filename = filename.strip()
if not filename.endswith('.csv'): filename = filename + '.csv'

# Folder name (from the basepath)
dirname = 'Input'
dirname = dirname.strip()

# File complete path
fpath = basepath / dirname 

# Delimiter used in the csv file
delimit = ','

# Initializing two 2d-arrays (number of time-steps during one day on axis = 0,
# number of typical days, i.e. months, on axis = 1).
n_months = 12

# Total time for each day and time-step (h)
time = 24
dt = 1
time_day = np.arange(0, time, dt)

n_timesteps = np.size(time_day)

# One will containing the sum of the production during a time-step of a certain months,
# over tha days of the months and the years considered, one will contain the total number
# of days summed in order to make an average
pv_production = np.zeros((n_timesteps, n_months))
pv_production_count_days = np.zeros((n_timesteps, n_months))


## Reading the original file

try:

    with open(fpath / filename, mode = 'r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter = delimit)

        # The row_before is a list containing the row that is before the current row
        # It is used to store the headers list
        row_before = []

        # The flag headers_flag is deactivated once that the header row has been read
        headers_flag = 1

        for row in csv_reader:

            # Different row are present in the beginning of the file that contains information
            # that are not needed , therefore rows which don't start with numerical values
            # or empty rows are skipped
            if not row == [] and row[0][0].isdigit():

                # When the first row with numerical values occurs, it means that the previous row was the
                # headers row, therefore it is stored in a list
                if headers_flag == 1:
                    headers = row_before
                    headers = [header.strip().lower().replace(' ', '_') for header in headers]
                    headers_flag = 0
                    
                # Each data-row is formatted as follows:
                # Column containing the date (header == 'time') aaaammdd:hhmm
                # Column containing the power (header == 'p'): ppp.pp
                # Values in the other columns are not needed
                # Notabene: the power is given in Watts, while the peak power is of 1 kWp,
                # since unit-production values are needed, the power is divided by 1000 Wp
                time_row = row[headers.index('time')]
                power = float(row[headers.index('p')])/1000

                # The values of the power are averaged for each month and for each time-step of the day,
                # among the days of each month and each year. In order to perform the average only the month and
                # the time-step (i.e. the hour) are needed. The months go from 1 (january) to 2 (december)
                # but the columns of the np.arrays go from 0 to 1, therefore 1 is subtracted to each month
                month = int(time_row[4:6]) - 1
                hour = int(time_row[9:11])
                
                # The total production for each time-step of each typical day is stored, and the number of days
                # considered too in order to perform an average afterwards
                pv_production[hour, month] += power
                pv_production_count_days[hour, month] += 1
              
            # Once that the headers have been stored, there's no more need to store the previous row
            if headers_flag == 1: row_before = row 

except: print('Unable to open this file')

# The pv_production array is given the proper shape: the time vector is added in the first column and the other columns
# are substituted by average values for the production
pv_production = np.column_stack((time_day, pv_production/pv_production_count_days))


## Storing the processed data in a .csv file
filename = 'pv_production_unit.csv'
fpath = basepath / dirname 
with open(fpath / filename , mode='w', newline='') as csv_file:
    csv_writer = csv.writer(csv_file, delimiter = ';', quotechar="'", quoting = csv.QUOTE_NONNUMERIC)

    csv_writer.writerow(['Time (h)'] + ['Month {} (kWh/h/kWp)'.format(i) for i in range(12)])

    for row in pv_production:
        csv_writer.writerow(row)