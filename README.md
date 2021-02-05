# EnergCommOptim

The repository contains a routine that optimizes the operation of a PV system with energy storage for fixed or variable (parametric) sizes for both of them, in the context of collective self-consumption and energy communities in Italy. PV production data are to be provided by the user (PVGIS database can be used), while consumption profiles are generated for an aggregate of households using probabilistic methods.

## Requirements

Codes included in this repository are written in Python 3, that is the only real requirement. They have been tested with Python 3.8 but also earlier version of Python 3 should work.
Python packages needed for running the methods are: pathlib, numpy, scipy, pulp, csv, tabulate, matplotlib.pyplot, math, random. All the other self-created methods used are present in this repository.

## Content of the repository

### Input/

In this folder, all the input _.csv_ files needed for the calculation are contained. Some of them must be updated from the user. Their name is properly formatted so that each methods knows which file to look at when certain data are needed. Particularly, the files are the followings. 


#### Optimization and energy assessment of the PV-storage system

These files should be updated from the user.

* *'pv_production_unit.csv'*: it contains the hourly unit production from the PV installation(s) in the given location, for typical days (one for each month of the year)

* *'battery_specs.csv'*: it contains the specifications for the battery (SOCmin, SOCmax, efficiencies,...)

* *'PVGIS_data.csv'*: it contains the hourly production from a PV installation in a given location for a number of years, as downloaded from PVGIS. It can be used to obtain *'pv_production_unit.csv'* if the latter is not provided from the user.


#### Generation of the load profiles

These files don't need to be updated from the user.

* *'eltdome_report.csv'*: it contains the attributes for all the appliances. It also contains, for each appliance, its "nickname", that is crucial for the correct loading of the other files.

* *'classenerg_report.csv'*: it contains the yearly energy consumption for each appliance, according to the different energy classes.

* *'coeff_matrix.csv'*: it contains the "user's behavior" coefficient for each appliance, according to the different seasons.

* Average daily load profile files for a type of appliance. They contain the time (in hours, from 00:00 to 23:59), generally with a resolution of 10 min, in the first column and the average power demand in Watt in the second column. Their name is formatted as follows: `'avg_loadprof' + '_' + app_nickname + '_' + seasons + '_' + day + '.csv'`. Where seasons indicates if the load profile is different according to the season (*'w'* stands for winter, *'s'* for summer, *'ap'* for autumn/spring) or the same for each season (*'sawp'*) and day indicates if the load profile is different according to the day in the week (*'wd'* stands for weekday and *'we'* for weekend day) or the same throughout the whole week (*'wde'*).

* Duty cycle files. They contain the time, with a resolution of 1 min, in the first columns and the power demand in Watt from the appliance in the second column. Their name is formatted as follows: `'dutycycle' + '_' + app_nickname + '.csv'`

#### Parameters/

This folder contains three .csv files where the parameters that are to be specified by the user for the simulation are stored.

* *'parameters.csv'*: it contains the user-defined values for some parameters that control the simulation.

#### Python files

* `aggregate_load_profiler_main.py`: this is the main. The load aggregated load profiles in the different seasons, and the energy consumptions over one year are computed and figures are generated. This is the only file that needs to be used directly by the user. The results are saved in Output/ (both files, in Files/, and figures in /Figures).

* `plot_generator.py`: this module contains all the methods for the creation of the figures showing the results.

* `parameters_input.py`: this module contains the method that is used in main for updating the paramters value to the user's keyboard input. This module is not directly used by the user. The parameters are saved in Parameters/.

* `house_load_profiler.py`: this module contains the method that computes the load profile for a household. This module can be used directly bu the user if the load profile for a single household is to be computed.

* `load_profiler.py`: this module contains the method that computes the load profile for a single appliance. This module can be used directly by the user if the load profile for a single appliance is to be computed.

* `cumulative_frequency.py`: this module contains a method to compute the cumulative frequency, starting from the frequency density. This module is not directely used by the user.

* `profile_interpolation.py`: this module contains the method that interpolates a given time-profile in order to change the time resolution. This module is not directely used by the user.

* `tictoc.py`: this module contains two methods that are a Python adaptation of Matalab's tic-toc functions.

* `levenshtein_distance.py`: this module contains a method that uses Levenshtein distance between two string to suggest the closest match to a word (user's input) and a list of words.

* `datareader.py`: this module contains different methods that properly read the various input files. This module is not directely used by the user.

#### Output/

This folder contains all the .csv files where the results from the various simulations are stored. The filename is formatted in order to give all the information needed about the location, season, day type of the simulation ,as well as the number of households considered and the energetic class of the appliances.
The folder also contains the subfolders Files/ and Figures/.
