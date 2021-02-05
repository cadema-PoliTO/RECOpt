# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 08:51:11 2020

@author: giamm
"""




import numpy as np
from tabulate import tabulate
from battery_optimization import battery_optimization
import matplotlib.pyplot as plt

def shared_energy_evaluator(time_dict, input_powers_dict, technologies_dict, auxiliary_dict, fixed_analysis_flag):



    ### Storing the given input in the proper variables

    ## Time discretization

    # # Total time of simulation (h) - for each typical day
    # time = time_dict['time']

    # Timestep for the simulation (h)
    dt = time_dict['dt']

    # # Vector of time, from 00:00 to 23:59, i.e. 24 h
    # time_sim = time_dict['time_sim']

    # Number of elements of the vector of time
    time_length = time_dict['time_length']


    ## Sizes and battery specficiations of the various technologies

    # # Grid maximum power (kW)
    # grid_power_max = technologies_dict['grid_power_max']

    # PV size (kW)
    pv_size = technologies_dict['pv_size']

    # # Battery size (kWh)
    # battery_size = technologies_dict['battery_size']

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
    n_days = auxiliary_dict['n_days']

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

    # In case of a fixed-size analysis, powers are stored for each typical day, in order to be plotted
    # This is done only if the fixed_analysis_flag is active in order to not slow down the procedure
    # storing such files if they are not needed
    if fixed_analysis_flag == 1:

        pv_production_month_day = np.zeros((time_length, n_months, n_days))
        grid_feed_month_day = np.zeros((time_length, n_months, n_days))
        grid_purchase_month_day = np.zeros((time_length, n_months, n_days))
        battery_charge_month_day = np.zeros((time_length, n_months, n_days))
        battery_discharge_month_day = np.zeros((time_length, n_months, n_days))
        battery_energy_month_day = np.zeros((time_length, n_months, n_days))
        shared_power_month_day = np.zeros((time_length, n_months, n_days))


    ## Performing the calculation for each month, for each day-type (i.e. for n_months*n_days = 24 typical days)

    # Since the keys in days (dict) are strings, this for loop is used to store their names (they are needed in the following)
    for day in days:
        iday = days[day][0]
        if iday == 0: weekday = day
        elif iday == 1: weekend = day

    # Running through the months
    for month in months:

        optimization_status[month] = {}

        # Unique number used to identify the month in the arrays elements
        mm = months[month]['id'][0]

        # In some cases the optimization method does not work properly and it is not able to perform the optimization.
        # When this happen, battery_optimization returns nans. In order not to discard a typical day in which the energy
        # values are nan from the total count, the energy is re-constructed using the adjacent typical day (same month)
        # or the adjacent month.
        # In more detail, if the energy values in a typical day are nan, an attempt is made to fill this values using
        # the energy values in the other typical day from the same month. Since weekdays are evaluated first, if energy
        # values for a weekday are nan, a flag is activated in order to fill them up at the next iteration (day type = weekend)
        weekday_nan_flag = 0

        # Running through the day-types (weekday or weekend day)
        for day in days:

            # Unique number used to identify the day-typein the arrays elements
            dd = days[day][0]

            # Initializing vectors where the powers are stored for the current month and day-type (i.e. typical day) for each time-step

            # Production from the PV (kWh/h)
            pv_production = np.zeros((time_length)) 

            # Consumption from the aggregate of households (kWh/h)
            consumption = np.zeros((time_length)) 

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
 
            # Evaluating total production, consumption for the typical day
            pv_production = pv_production_month[:, mm]  
            consumption = consumption_month_day[:, mm, dd]
            
            # Using the method battery_optimization from the module battery_optimization.py to run the MILP
            # optimization in order to determine the battery usage strategy in the typical day. Powers and 
            # energy stored in the battery at each time-step are returned.

            # optimization_status[month][day], \
            # grid_feed, grid_purchase, battery_charge, battery_discharge, battery_energy = \
            # battery_optimization(pv_production, pv_available, net_load, time_dict, technologies_dict)

            optimization_status[month][day], \
            shared_power, grid_feed, grid_purchase, battery_charge, battery_discharge, battery_energy = \
            battery_optimization(pv_production, consumption, time_dict, technologies_dict)
 
            # # Uncomment to check that the constraints and equations defined in the problem are actually respected
            # tol = 1e-4
            # if (np.any(net_load + grid_feed + battery_charge - pv_available - grid_purchase - battery_discharge > tol)): print('ops, node')
            # if (np.any(grid_feed > grid_power_max + tol) or np.any(grid_feed > grid_power_max + tol)): print('ops, grid power')
            # if (np.any(grid_feed > pv_available + tol)): print('ops, grid feed')
            # if (np.any(battery_charge > pv_available + tol)): print('ops, battery_charge'); print(battery_charge - pv_available)
        
            # Evaluating the instanteneous shared energy, according to the definition in the Decree Law 162/2019
            # shared_power = np.minimum((pv_production + battery_discharge), (consumption + battery_charge))

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


            # In case nan are returned by battery_optimization, the contribution of the current typical day to the energy
            # in the current month is null (due to np.nansum), therefore the value is to be fixed
            if np.any(np.isnan(shared_power)): 

                # print('{}, {} energy is nan'.format(month, day))
                # print('Shared energy before fixing: {}'.format(shared_energy[mm]))

                # If the day is weekday (that is evaluated first, for each month) it is needed to wait unitl the next
                # iteration, therefore a flag is activated
                if day == weekday:
                    weekday_nan_flag = 1

                # If the day is weekend day, the values can be fixed using the values from the weekday of the same month
                elif day == weekend and weekday_nan_flag == 0:
                    n_days_weekday = months[month]['days_distr'][weekday]
                    grid_feed_energy[mm] += grid_feed_energy[mm]/n_days_weekday*number_of_days
                    grid_purchase_energy[mm] += grid_purchase_energy[mm]/n_days_weekday*number_of_days
                    battery_charge_energy[mm] += battery_charge_energy[mm]/n_days_weekday*number_of_days
                    battery_discharge_energy[mm] += battery_discharge_energy[mm]/n_days_weekday*number_of_days
                    shared_energy[mm] += shared_energy[mm]/n_days_weekday*number_of_days

                    # print('Shared energy after fixing: {}'.format(shared_energy[mm]))

                # If both the weekday and weekend day have retuned nan, the energy for the month is zero. Anyway
                # it is set to nan in order to fix it at the end of the for loops, interpolating between adjacent months
                elif day == weekend and weekday_nan_flag == 1:
                    grid_feed_energy[mm] = np.nan
                    grid_purchase_energy[mm] = np.nan
                    battery_charge_energy[mm] = np.nan
                    battery_discharge_energy[mm] = np.nan
                    shared_energy[mm] = np.nan

                    # print('Shared energy after fixing: {}'.format(shared_energy[mm]))
        
            # If the day is weekend and the flag from the previous weekday (same month) has been activated,
            # the energy is fixed
            elif day == weekend and weekday_nan_flag == 1 :

                # print('Shared energy before fixing: {}'.format(shared_energy[mm]))

                n_days_weekday = months[month]['days_distr'][weekday]
                grid_feed_energy[mm] += grid_feed_energy[mm]/number_of_days*n_days_weekday
                grid_purchase_energy[mm] += grid_purchase_energy[mm]/number_of_days*n_days_weekday
                battery_charge_energy[mm] += battery_charge_energy[mm]/number_of_days*n_days_weekday
                battery_discharge_energy[mm] += battery_discharge_energy[mm]/number_of_days*n_days_weekday
                shared_energy[mm] += shared_energy[mm]/number_of_days*n_days_weekday

                # print('Shared energy after fixing: {}'.format(shared_energy[mm]))

            # In case of a fixed-size analysis, powers are stored for each typical day, in order to be plotted
            # This is done only if the fixed_analysis_flag is active in order to not slow down the procedure
            # storing such files if they are not needed
            # Nota bene: fixed-size analysis means that both the pv and battery size are fixed (to one value)
            # it has nothing to do with the fixing procedure operated for nan values
            if fixed_analysis_flag == 1:

                pv_production_month_day[:, mm, dd] = pv_production
                grid_feed_month_day[:, mm, dd] = grid_feed
                grid_purchase_month_day[:, mm, dd] = grid_purchase
                battery_charge_month_day[:, mm, dd] = battery_charge
                battery_discharge_month_day[:, mm, dd] = battery_discharge
                battery_energy_month_day[:, mm, dd] = battery_energy
                shared_power_month_day[:, mm, dd] = shared_power


    # If both the weekday and weekend returned nan values, the energy for the month has been set to nan
    # Now it is fixed using the energy values from the adjacent months (interpolating the values)
    if np.any(np.isnan(shared_energy)):

        nans, fnan = np.isnan(grid_feed_energy), lambda z: z.nonzero()[0]
        grid_feed_energy[nans]= np.interp(fnan(nans), fnan(~nans), grid_feed_energy[~nans], period = n_months)

        nans, fnan = np.isnan(grid_purchase_energy), lambda z: z.nonzero()[0]
        grid_purchase_energy[nans]= np.interp(fnan(nans), fnan(~nans), grid_purchase_energy[~nans], period = n_months)

        nans, fnan = np.isnan(battery_charge_energy), lambda z: z.nonzero()[0]
        battery_charge_energy[nans]= np.interp(fnan(nans), fnan(~nans), battery_charge_energy[~nans], period = n_months)

        nans, fnan = np.isnan(battery_discharge_energy), lambda z: z.nonzero()[0]
        battery_discharge_energy[nans]= np.interp(fnan(nans), fnan(~nans), battery_discharge_energy[~nans], period = n_months)

        nans, fnan = np.isnan(shared_energy), lambda z: z.nonzero()[0]
        shared_energy[nans]= np.interp(fnan(nans), fnan(~nans), shared_energy[~nans], period = n_months)


    results = {
        'optimization_status': optimization_status,
        'pv_production_energy': pv_production_energy,
        'consumption_energy': consumption_energy,
        'grid_feed_energy': grid_feed_energy,
        'grid_purchase_energy': grid_purchase_energy,
        'battery_charge_energy': battery_charge_energy,
        'battery_discharge_energy': battery_discharge_energy,
        'shared_energy': shared_energy,       
        }

    # In case of a fixed analysis, powers are stored for each typical day, in order to be plotted
    # This is done only if the fixed_analysis_flag is active in order to not slow down the procedure
    # storing such files if they are not needed
    if fixed_analysis_flag == 1:
        results['pv_production_month_day'] = pv_production_month_day
        results['consumption_month_day'] = consumption_month_day
        results['grid_feed_month_day'] = grid_feed_month_day
        results['grid_purchase_month_day'] = grid_purchase_month_day
        results['battery_charge_month_day'] = battery_charge_month_day
        results['battery_discharge_month_day'] = battery_discharge_month_day
        results['battery_energy_month_day'] = battery_energy_month_day
        results['shared_power_month_day'] = shared_power_month_day

    return results




        
    





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