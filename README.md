# RECOpt

## Energetic evaluation and optimization of renewable energy communities

The repository contains a tool that optimizes the operation of a PV system with energy storage for fixed or variable (parametric) sizes for both of them, in the context of collective self-consumption and energy communities in Italy. PV production data are to be provided by the user (PVGIS database can be used), while consumption profiles are generated for an aggregate of households using probabilistic methods.

## Requirements

Codes included in this repository are written in Python 3, that is the only real requirement. They have been tested with Python 3.8 but also earlier version of Python 3 should work.
Python packages needed for running the methods are: pathlib, numpy, scipy, pulp, csv, tabulate, matplotlib.pyplot, math, random. All the other self-created methods needed for the tool to work are contained in this repository.

## Content of the repository

### Input/

In this folder, all the input _.csv_ files needed for the calculation are contained. Some of them must be updated from the user. Their name is properly formatted so that each methods knows which file to look at when certain data are needed. Particularly, the files are the followings. 


#### Optimization and energy assessment of the PV-storage system

These files should be updated by the user.

* *'pv_production_unit.csv'*: it contains the hourly unit production from the PV installation(s) in the given location, for typical days (one for each month of the year)

* *'battery_specs.csv'*: it contains the specifications for the battery (SOCmin, SOCmax, efficiencies,...)

* *'PVGIS_data.csv'*: it contains the hourly production from a PV installation in a given location for a number of years, as downloaded from PVGIS. It can be used to obtain *'pv_production_unit.csv'* if the latter is not provided from the user.


#### Generation of the load profiles

These files don't need to be updated by the user.

* *'eltdome_report.csv'*: it contains the attributes for all the appliances. It also contains, for each appliance, its "nickname", that is crucial for the correct loading of the other files.

* *'classenerg_report.csv'*: it contains the yearly energy consumption for each appliance, according to the different energy classes.

* *'coeff_matrix.csv'*: it contains the "user's behavior" coefficient for each appliance, according to the different seasons.

* Average daily load profile files for a type of appliance. They contain the time (in hours, from 00:00 to 23:59), generally with a resolution of 10 min, in the first column and the average power demand in Watt in the second column. Their name is formatted as follows: `'avg_loadprof' + '_' + app_nickname + '_' + seasons + '_' + day + '.csv'`. Where seasons indicates if the load profile is different according to the season (*'w'* stands for winter, *'s'* for summer, *'ap'* for autumn/spring) or the same for each season (*'sawp'*) and day indicates if the load profile is different according to the day in the week (*'wd'* stands for weekday and *'we'* for weekend day) or the same throughout the whole week (*'wde'*).

* Duty cycle files. They contain the time, with a resolution of 1 min, in the first columns and the power demand in Watt from the appliance in the second column. Their name is formatted as follows: `'dutycycle' + '_' + app_nickname + '.csv'`


### Python files

* `main.py`: this is the only file that should be directly used by the user. No manual modifications should be made, e.g. to change some parameters; the user just needs to make it run. Parameters are updated from keyboard and stored in Parameters/. Results, both _.csv_ files and _.png_ figures, are stored in Outputs/ for each in simulation, respectively, in Files/ and Figures/.

* `pvgis_to_csv.py`: the module processes the data about hourly production from the PV, contained in a _.csv_ file downloaded from PVGIS. This module should be used from the user if the *'pv_production_unit.csv'* file is not directly provided.

* `shared_energy_evaluator.py`: the module contains a method that evaluates the performance (shared energy, and other quantities) for a given configuration (number of households, size of the PV and battery systems) in one year, using a number of typical days (two for each month, both week-day and weekend-day). This module is not directly used by the user. 

* `battery_optimization.py`: the module contains a method that optimizes the operation of the battery in one day, once that the production from the pv and the consumption from the households are given. At the moment, if the user wants to change the objective of the optimization, this should be done here, manually. This module is not directly used by the user. 

* `aggregate_load_profiler.py`: the module contains a method that generates the aggregated load profiles in different typical days (two for each season, both week-day and weekend-day) for a number of households. If the user chooses so when running the simulation, detailed files and figures about the load profiles and the energy consumption from the electric appliances are generated and saved in Output/. If the user wants to change some very specific parameters about the generation of the load profiles, this should be done here. Normally the user does not need to use this module, but it can be used to test the generation of the load profile for the aggregate of household.

* `house_load_profiler.py`: this module contains a method that computes the load profile for a household in a given typical day. Normally the user does not need to use this module, but it can be used to test the generation of the load profile for a single household.

* `load_profiler.py`: this module contains a method that computes the load profile for a single appliance in a given typical day. Normally the user does not need to use this module, but it can be used to test the generation of the load profile for a single appliance.

* `load_profile_aggregator_trapz.py`: this module contains a method that aggregates some profiles using a different time-step. This module is not directly used by the user. 

* `plot_generator.py`: this module contains all the methods for the creation of the figures showing the results. This module is not directly used by the user. 

* `parameters_input.py`: this module contains the methods that are used in main for updating the paramters value to the user's keyboard input. The parameters are saved in Parameters/. This module is not directly used by the user. 

* `levenshtein_distance.py`: this module contains a method that uses Levenshtein distance between two string to suggest the closest match to a word (user's input) and a list of words. This module is not directly used by the user. 

* `datareader.py`: this module contains different methods that properly read the various input files. This module is not directely used by the user.

* `tictoc.py`: this module contains two methods that are a Python adaptation of Matalab's tic-toc functions. This module is not directly used by the user. 




