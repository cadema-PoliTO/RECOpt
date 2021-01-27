# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 08:51:11 2020

@author: giamm
"""




import numpy as np
from tabulate import tabulate
from battery_optimization import battery_optimization
import matplotlib.pyplot as plt

def shared_energy_evaluator(time_dict, input_powers_dict, technologies_dict, auxiliary_dict):



    ### Storing the given input in the proper variables

    ## Time discretization

    # # Total time of simulation (h) - for each typical day
    # time = time_dict['time']

    # Timestep for the simulation (h)
    dt = time_dict['dt']

    # Vector of time, from 00:00 to 23:59, i.e. 24 h
    time_sim = time_dict['time_sim']

    # Number of elements of the vector of time
    time_length = time_dict['time_length']


    ## Sizes and battery specficiations of the various technologies

    # # Grid maximum power (kW)
    # grid_power_max = technologies_dict['grid_power_max']

    # PV size (kW)
    pv_size = technologies_dict['pv_size']

    # Battery size (kW)
    battery_size = technologies_dict['battery_size']

    # # Battery specifications 
    # battery_specs = technologies_dict['battery_specs']


    ## Input powers for each month/month and day

    # PV production during each month (kWh/h/kWp)
    pv_production_unit = input_powers_dict['pv_production_unit']

    # Actual production from the PV (kWh/h)
    pv_production_month = pv_production_unit*pv_size

    # Consumption from the households during each month and day-type (kWh/h)
    consumption_month_day = input_powers_dict['consumption_month_day']


    ## Auxiliary variable (seasons/month/days/days_distr dictionaries)

    # # Seasons dict (how many seasons)
    # seasons = auxiliary_dict['seasons']
    # n_seasons = auxiliary_dict['n_seasons']

    # Months dict (how many months and which season do they belong to)
    months = auxiliary_dict['months']
    n_months = auxiliary_dict['n_months']

    # Days dict (how many day-types: weekday and weekend day)
    days = auxiliary_dict['days']
    # n_days = auxiliary_dict['n_days']

    # Days distribution (how many days for each day-type during each month)
    # days_distr = auxiliary_dict['days_distr']



    ### Starting the evaluation

    # Optimization status for each month and day-type

    optimization_status = {}

    # Monthly values for the energy are stored 

    # Energy production from the PV (kWh/month)
    pv_production_energy = np.zeros((n_months))

    # Energy consumption from the aggregate of households (kWh/month)
    consumption_energy = np.zeros((n_months))

    # Excess energy from the PV that is fed into the grid (kWh/month)
    # Nota bene: since a virtual scheme is considered, all the energy that is produced is fed into the grid
    # Therefore it is important to keep in mind the difference between the energy that is fed into the grid 
    # to be consumed locally and the one that is fed because there is an excess. The term "feed" will be used 
    # for the latter only
    grid_feed_energy = np.zeros((n_months))

    # Consumption that is satisfied purchasing energy from the grid (kWh/month)
    # Nota bene: similarly as above, all the consumption is satisfied taking energy from the grid.
    # Anyway, in order to differentiate the energy that is shared from the deficit energy that is
    # actually satisfied from the grid, the term "purchase" is used.
    grid_purchase_energy = np.zeros((n_months))

    # Excess energy injected into the battery (kWh/month)
    battery_charge_energy = np.zeros((n_months))

    # Consumption that is satisfied taking energy from the battery (kWh/month)
    battery_discharge_energy = np.zeros((n_months))

    # Shared energy (according to the definition provided by the Decree Law 162/2019 - Milleproroghe) (kWh/month)
    shared_energy = np.zeros((n_months))


    ## Performing the calculation for each month, for each day-type (i.e. for n_months*n_days = 24 typical days)

    # Running through the months
    for month in months:

        optimization_status[month] = {}

        # Unique number used to identify the month in the arrays elements
        mm = months[month]['id'][0]

        # # Uncomment to plot the power values for this typical day
        # plot_params = {
        # 'time_scale': 'h',
        # 'power_scale': 'kW',
        # 'energy_scale': 'MWh',
        # 'figsize': (297/25.4 , 420/25.4),
        # 'orientation': 'horizontal',
        # 'font_small': 14,
        # 'font_medium': 16,
        # 'font_large': 18,
        # }

        # # Figure setup: figure size and orientation, font-sizes 
        # figsize = plot_params['figsize']
        # orientation = plot_params['orientation']

        # if orientation == 'horizontal': figsize = figsize[::-1]

        # fontsize_title = plot_params['font_large']
        # fontsize_legend = plot_params['font_medium']
        # fontsize_labels = plot_params['font_medium']
        # fontsize_text = plot_params['font_medium']
        # fontsize_ticks = plot_params['font_small']
        # fontsize_pielabels = plot_params['font_small']

        # # Creating a figure with multiple subplots, with two rows (one for each type of day)
        # fig, ax = plt.subplots(2, 1, sharex = False, sharey = False, figsize = figsize)
        
        # # suptitle = 'Results in {}'.format(month)
        # # fig.suptitle(suptitle, fontsize = fontsize_title, fontweight = 'bold')
        # fig.subplots_adjust(left = 0.1, bottom = 0.1, right = 0.9, top = 0.85, wspace = None, hspace = 0.3)

        # Running through the day-types (weekday or weekend day)
        for day in days:

            # Unique number used to identify the day-typein the arrays elements
            dd = days[day][0]

            # Initializing vectors where the powers are stored for the current month and day-type (i.e. typical day) for each time-step

            # Production from the PV (kWh/h)
            pv_production = np.zeros((time_length)) 

            # Consumption from the aggregate of households (kWh/h)
            consumption = np.zeros((time_length)) 

            # Net load, i.e. consumption that cannot be satisfied by the PV production (kWh/h)
            net_load = np.zeros((time_length))  

            # Available production, i.e. production that remains after satisfying the local demand (kWh/h)
            pv_available = np.zeros((time_length)) 

            # Excess power (available) that charges the battery (kWh/h)
            battery_charge= np.zeros((time_length)) 

            # Power that discharges the battery (kWh/h)
            battery_discharge = np.zeros((time_length)) 

            # Excess power that is fed into the grid (kWh/h)
            grid_feed = np.zeros((time_length)) 

            # Deficit power that is taken from the grid (kWh/h)
            grid_purchase = np.zeros((time_length)) 

            # Instanteously shared energy (kWh/h)
            shared_power = np.zeros((time_length)) 

            # Energy that is stored in the energy at each time-step (kWh)
            battery_energy = np.zeros((time_length)) 
 

            # Evaluating production, consumption and net production/load in the current typical day
            pv_production = pv_production_month[:, mm]  
            consumption = consumption_month_day[:, mm, dd]

            pv_available = pv_production - consumption
            # net_load[pv_available < 0] = -pv_available[pv_available < 0]
            pv_available[pv_available < 0]= 0
            net_load = consumption - pv_production
            
            # Using the method battery_optimization from the module battery_optimization.py to run the MILP
            # optimization in order to determine the battery usage strategy in the typical day. Powers and 
            # energy stored in the battery at each time-step are returned.

            optimization_status[month][day], \
            grid_feed, grid_purchase, battery_charge, battery_discharge, battery_energy = \
            battery_optimization(pv_available, net_load, time_dict, technologies_dict)

            # from pyxems13 import compute_cons
            # nettt_load = consumption - pv_production
            # battery_charge, battery_discharge, grid_feed, grid_purchase, battery_energy = \
            # compute_cons(pv_available, nettt_load, battery_size, time_length - 1, time_sim + 1, dt)
            # battery_charge, battery_discharge, grid_feed, grid_purchase, battery_energy = \
            # np.asarray(battery_charge), np.asarray(battery_discharge), np.asarray(grid_feed), np.asarray(grid_purchase), np.asarray(battery_energy)
 

            # # Uncomment to check that the constraints and equations defined in the problem are actually respected
            # tol = 1e-4
            # if (np.any(net_load + grid_feed + battery_charge - pv_available - grid_purchase - battery_discharge > tol)): print('ops, node')
            # if (np.any(grid_feed > grid_power_max + tol) or np.any(grid_feed > grid_power_max + tol)): print('ops, grid power')
            # if (np.any(grid_feed > pv_available + tol)): print('ops, grid feed')
            # if (np.any(battery_charge > pv_available + tol)): print('ops, battery_charge'); print(battery_charge - pv_available)
        
            # Evaluating the instanteneous shared energy, according to the definition in the Decree Law 162/2019
            shared_power = np.minimum((pv_production + battery_discharge), (consumption + battery_charge))

            # Storing the energies for this typical day in the monthly energies, this is done by summing the power during each time-step
            # (and multipling by the time-step in order to get the energy) and multiplying by the total number of days of this type
            number_of_days = months[month]['days_distr'][day]

            pv_production_energy[mm] += np.nansum(pv_production)*dt*number_of_days
            grid_feed_energy[mm] += np.nansum(grid_feed)*dt*number_of_days
            grid_purchase_energy[mm] += np.nansum(grid_purchase)*dt*number_of_days
            consumption_energy[mm] += np.nansum(consumption)*dt*number_of_days
            battery_charge_energy[mm] += np.nansum(battery_charge)*dt*number_of_days
            battery_discharge_energy[mm] += np.nansum(battery_discharge)*dt*number_of_days
            shared_energy[mm] += np.nansum(shared_power)*dt*number_of_days

            # # Uncomment to plot the power values for this typical day

            # # Plotting the value for the typical day
            # ax[dd].plot(time_sim + dt/2, pv_production, label = 'pv_production')
            # ax[dd].plot(time_sim + dt/2, consumption, label = 'consumption')
            # ax[dd].plot(time_sim + dt/2, grid_feed, label = 'grid_feed')
            # ax[dd].plot(time_sim + dt/2, grid_purchase, label = 'grid_purchase')
            # ax[dd].plot(time_sim + dt/2, battery_charge, label = 'battery_charge')
            # ax[dd].plot(time_sim + dt/2, battery_discharge, label = 'battery_discharge')

            # ax[dd].bar(time_sim, shared_power, color = 'k', width = dt, align = 'edge', label = 'shared power', alpha = 0.3)
            
            # title = '{}, {}'.format(month.capitalize(), day)
            # ax[dd].set_title(title, fontsize = fontsize_title)

            # axtw = ax[dd].twinx()
            # axtw.plot(time_sim + dt/2, battery_energy/battery_size*100, 'r--', label = 'energy battery')
            
            # # Making the figure look properly
       
            # ax[dd].set_xlabel('Time ({})'.format('h'), fontsize = fontsize_labels)
            # ax[dd].set_ylabel('Power ({})'.format('kW'), fontsize = fontsize_labels)
            # ax[dd].set_xlim([time_sim[0], time_sim[-1]])
            # # Set one tick each hour on the x-axis
            # ax[dd].set_xticks(list(time_sim[: : int(dt)]))
            # ax[dd].tick_params(axis ='both', labelsize = fontsize_ticks)
            # ax[dd].tick_params(axis ='x', labelrotation = 0)
            # ax[dd].grid()
            # ax[dd].legend(loc = 'upper left', fontsize = fontsize_legend, ncol = 2)

            # axtw.set_ylabel('SOC ({})'.format('%'), fontsize = fontsize_labels)
            # axtw.legend(loc = 'upper right', fontsize = fontsize_legend, ncol = 2)
            # axtw.set_ylim([0, 100])

  

 




    # index_self_suff_month = shared_energy/(consumption_energy + battery_charge_energy)*100
    # index_self_cons_month = shared_energy/(pv_production_energy + battery_discharge_energy)*100

    # index_self_suff_year =  np.sum(shared_energy)/np.sum(consumption_energy + battery_charge_energy)*100
    # index_self_cons_year = np.sum(shared_energy)/np.sum(pv_production_energy + battery_discharge_energy)*100


    results = {
        'optimization_status': optimization_status,
        'pv_production_energy': pv_production_energy,
        'consumption_energy': consumption_energy,
        'grid_feed_energy': grid_feed_energy,
        'grid_purchase_energy': grid_purchase_energy,
        'battery_charge_energy': battery_charge_energy,
        'battery_discharge_energy': battery_discharge_energy,
        'shared_energy': shared_energy,
        # 'index_self_suff_month': index_self_suff_month,
        # 'index_self_cons_month': index_self_cons_month,
        # 'index_self_suff_year': index_self_suff_year,
        # 'index_self_cons_year': index_self_cons_year,        
        }

    return results



   



    # [print('{0:s}: ISS: {1:.2f}%, ISC: {2:.2f}%'.format(month.capitalize(), index_self_suff_month[months[month]['id'][0]], index_self_cons_month[months[month]['id'][0]])) for month in months]

    # print('\nYear: ISS: {0:.2f}%, ISC: {1:.2f}%'.format(index_self_suff_year, index_self_cons_year))


        
    





    # ### Starting the evaluation

    # # Monthly values for the energy are stored in propero variables, divided by day-type (weekday or weekend day)

    # # Energy production from the PV (kWh/month/day-type)
    # pv_production_energy = np.zeros((n_months, n_days))

    # # Energy consumption from the aggregate of households (kWh/month/day-type)
    # consumption_energy = np.zeros((n_months, n_days))

    # # Excess energy from the PV that is fed into the grid and sold (kWh/month/day-type)
    # # Nota bene: since a virtual scheme is considered, all the energy that is produced is fed into the grid
    # # Therefore it is important to keep in mind the difference between the energy that is fed into the grid 
    # # to be consumed locally and the one that is fed because there is an excess. The term "feed" will be used 
    # # for the latter only
    # grid_feed_energy = np.zeros((n_months, n_days))

    # # Consumption that is satisfied purchasing energy from the grid (kWh/month/day-type)
    # # Nota bene: similarly as above, all the consumption is satisfied taking energy from the grid.
    # # Anyway, in order to differentiate the energy that is shared from the deficit energy that is
    # # actually satisfied from the grid, the term "purchase" is used.
    # grid_purchase_energy = np.zeros((n_months, n_days))

    # # Excess energy injected into the battery (kWh/month/day-type)
    # battery_charge_energy = np.zeros((n_months, n_days))

    # # Consumption that is satisfied taking energy from the battery (kWh/month/day-type)
    # battery_discharge_energy = np.zeros((n_months, n_days))

    # # Shared energy (according to the definition provided by the Decree Law 162/2019 - Milleproroghe) (kWh/month/day-type)
    # shared_energy = np.zeros((n_months, n_days))

    # # The calculation is performed for each month, for each day-type (i.e. for n_months*n_days = 24 typical days)
    # for month in months:

    #     # Unique number used to identify the month in the arrays elements
    #     mm = months[month]['id'][0]

    #     # Initializing vectors where the powers are stored for the current month, for both day-types, for each time-step

    #     # Production from the PV (kWh/h)
    #     pv_production = np.zeros((time_length, n_days)) 

    #     # Consumption from the aggregate of households (kWh/h)
    #     consumption = np.zeros((time_length, n_days))

    #     # Net load, i.e. consumption that cannot be satisfied by the PV production (kWh/h)
    #     net_load = np.zeros((time_length, n_days)) 

    #     # Available production, i.e. production that remains after satisfying the local demand (kWh/h)
    #     pv_available = np.zeros((time_length, n_days))

    #     # Excess power (available) that charges the battery (kWh/h)
    #     battery_charge= np.zeros((time_length, n_days))

    #     # Power that discharges the battery (kWh/h)
    #     battery_discharge = np.zeros((time_length, n_days))

    #     # Excess power that is fed into the grid (kWh/h)
    #     grid_feed = np.zeros((time_length, n_days))

    #     # Deficit power that is taken from the grid (kWh/h)
    #     grid_purchase = np.zeros((time_length, n_days))

    #     # Instanteously shared energy (kWh/h)
    #     shared_power = np.zeros((time_length, n_days))

    #     # Energy that is stored in the energy at each time-step (kWh)
    #     battery_energy = np.zeros((time_length, n_days))

    #     # Creating a figure with multiple subplots, with two rows (one for each type of day)
    #     fig, ax = plt.subplots(2, 1, sharex = False, sharey = False, figsize = figsize)
        
    #     # suptitle = 'Results in {}'.format(month)
    #     # fig.suptitle(suptitle, fontsize = fontsize_title, fontweight = 'bold')
    #     fig.subplots_adjust(left = 0.1, bottom = 0.1, right = 0.9, top = 0.85, wspace = None, hspace = 0.3)

        
    #     for day in days:

    #         dd = days[day][0]


    #         pv_production[:, dd]  = pv_production_month[:, mm]  
    #         consumption[:, dd] = consumption_month_day[:, mm, dd]

    #         pv_available[:, dd] = pv_production[:, dd] - consumption[:, dd]
    #         net_load[pv_available[:, dd] < 0, dd] = -pv_available[pv_available[:, dd] < 0, dd]
    #         pv_available[pv_available[:, dd] < 0, dd]= 0
            
    #         grid_feed[:, dd], grid_purchase[:, dd], battery_charge[:, dd], battery_discharge[:, dd], battery_energy[:, dd] = \
    #         battery_optimization(pv_available[:, dd], net_load[:, dd], time_dict, technologies_dict)
        
    #         shared_power[:, dd] = np.minimum((pv_production[:, dd] + battery_discharge[:, dd]), (consumption[:, dd] + battery_charge[:, dd]))


    #         tol = 1e-4
    #         if (np.any(net_load[:, dd] + grid_feed[:, dd] + battery_charge[:, dd] - pv_available[:, dd] - grid_purchase[:, dd] - battery_discharge[:, dd] > tol)): print('ops, node')
    #         if (np.any(grid_feed[:, dd] > grid_power_max + tol) or np.any(grid_feed[:, dd] > grid_power_max + tol)): print('ops, grid power')
    #         if (np.any(grid_feed[:, dd] > pv_available[:, dd] + tol)): print('ops, grid feed')
    #         if (np.any(battery_charge[:, dd] > pv_available[:, dd] + tol)): print('ops, battery_charge')

    #         number_of_days = months[month]['days_distr'][day]

    #         pv_production_energy[mm][dd] = np.sum(pv_production)*dt*number_of_days
    #         grid_feed_energy[mm][dd] = np.sum(grid_feed)*dt*number_of_days
    #         grid_purchase_energy[mm][dd] = np.sum(grid_purchase)*dt*number_of_days
    #         consumption_energy[mm][dd] = np.sum(consumption)*dt*number_of_days
    #         battery_charge_energy[mm][dd] = np.sum(battery_charge)*dt*number_of_days
    #         battery_discharge_energy[mm][dd] = np.sum(battery_discharge)*dt*number_of_days
    #         shared_energy[mm][dd] = np.sum(shared_power)*dt*number_of_days


    #         ax[dd].plot(time_sim + dt/2, pv_available[:, dd], label = 'pv_avaialble')
    #         ax[dd].plot(time_sim + dt/2, net_load[:, dd], label = 'net_load')
    #         ax[dd].plot(time_sim + dt/2, grid_feed[:, dd], label = 'grid_feed')
    #         ax[dd].plot(time_sim + dt/2, grid_purchase[:, dd], label = 'grid_purchase')
    #         ax[dd].plot(time_sim + dt/2, battery_charge[:, dd], label = 'battery_charge')
    #         ax[dd].plot(time_sim + dt/2, battery_discharge[:, dd], label = 'battery_discharge')

    #         ax[dd].bar(time_sim, shared_power[:, dd], color = 'k', width = dt, align = 'edge', label = 'shared power', alpha = 0.3)
            
    #         title = '{}, {}'.format(month.capitalize(), day)
    #         ax[dd].set_title(title, fontsize = fontsize_title)

    #         axtw = ax[dd].twinx()
    #         axtw.bar(time_sim, battery_energy[:, dd]/battery_capacity*100, color = 'y', width = dt, align = 'edge', label = 'shared power', alpha = 0.3)
            

    #     pv_production_energy_month = np.sum(pv_production_energy, axis = 1)
    #     grid_feed_energy_month = np.sum(grid_feed_energy, axis = 1)
    #     grid_purchase_energy_month = np.sum(grid_purchase_energy, axis = 1)
    #     consumption_energy_month= np.sum(consumption_energy, axis = 1)
    #     battery_charge_energy_month = np.sum(battery_charge_energy, axis = 1)
    #     battery_discharge_energy_month = np.sum(battery_discharge_energy, axis = 1)
    #     shared_energy_month = np.sum(shared_energy, axis = 1)

    #     index_self_suff_month = shared_energy_month/consumption_energy_month*100
    #     index_self_cons_month = shared_energy_month/pv_production_energy_month*100

    #     index_self_suff_year =  np.sum(shared_energy_month)/np.sum(consumption_energy_month)*100
    #     index_self_cons_year = np.sum(shared_energy_month)/np.sum(pv_energy_month)*100

    #     [print('{0:s}: ISS: {1:.2f}%, ISC: {2:.2f}%'.format(month.capitalize(), index_self_suff_month[months[month]['id'][0]], index_self_cons_month[months[month]['id'][0]])) for month in months]

    #     print('\nYear: ISS: {0:.2f}%, ISC: {1:.2f}%'.format(index_self_suff_year, index_self_cons_year))
        
        
        
    #     ##
    #     # Making the figure look properly
    #     for axi in ax.flatten():
    #         axi.set_xlabel('Time ({})'.format('h'), fontsize = fontsize_labels)
    #         axi.set_ylabel('Power ({})'.format('kW'), fontsize = fontsize_labels)
    #         axi.set_xlim([time_sim[0], time_sim[-1]])
    #         # Set one tick each hour on the x-axis
    #         axi.set_xticks(list(time_sim[: : int(dt)]))
    #         axi.tick_params(axis ='both', labelsize = fontsize_ticks)
    #         axi.tick_params(axis ='x', labelrotation = 0)
    #         axi.grid()
    #         axi.legend(loc = 'upper left', fontsize = fontsize_legend, ncol = 2)
        

    # plt.show()