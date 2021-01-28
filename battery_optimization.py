import numpy as np
from pulp import *


###############################################################################


# This module contains the method that performs the optimization problem defining
# the operational strategy of the battery


###############################################################################

def battery_optimization(power_available, net_load, time_dict, technologies_dict):
    '''
    The method optimizes the operation of the battery during one day.
    Input:
        power_available: 1d-array, containing the value of the excess power at each time-step
        net_load: 1d-array, containing the value of the net load (consumption - production) at each time-step
        time_dict: dict, containing all the elements needed for the time discretization
        technologies_dict: dict, containing all the information about the technologies involved (PV, battery, grid)

    Output:
        optimization_status: str, showing the status of the optimization
        grid_feed: 1d-array, containing the excess power fed into the grid at each timestep
        grid_purchase: 1d-array, containing the deficit power purchased the grid at each timestep
        battery_charge: 1d-array, containing the excess power used to charge the battery at each timestep
        battery_discharge: 1d-array, containing the deficit power taken from the battery at each timestep
        battery_energy: 1d-array, containing the amount of power stored in the battery 
    '''



    ### Storing the given input in the proper variables

    ## Components activation
    # It is possible to exclude some components from the simulation

    grid_purchase_flag = 1
    grid_feed_flag = 1
    battery_flag = 1


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

    # Grid maximum power (kW)
    grid_power_max = technologies_dict['grid_power_max']

    # # PV size (kW)
    # pv_size = technologies_dict['pv_size']

    # Battery size/capacity (kWh)
    battery_size = technologies_dict['battery_size']
    battery_capacity = battery_size

    # Battery specifications 
    battery_specs = technologies_dict['battery_specs']
     
    
    ## Battery specifications
    
    # Maximum and minimum states of charge (SOC) (%)
    SOCmax = battery_specs['SOC_max']
    SOCmin = battery_specs['SOC_min']

    # Minimum time of charge/discharge (on which the maximum discharge power depends) (h)
    t_cd_min = battery_specs['t_cd_min']

    # Charge, discharge and self-discharge efficiencies (-)
    eta_charge = battery_specs['eta_charge']
    eta_discharge = battery_specs['eta_discharge']
    eta_self_discharge = battery_specs['eta_self_discharge']

    # Maximum and minimum energy in the battery (kWh)
    battery_energy_max = SOCmax*battery_capacity 
    battery_energy_min = SOCmin*battery_capacity 

    # Maximum power of discharge and charge (kW)
    battery_discharge_pmax = battery_capacity*(SOCmax-SOCmin)/t_cd_min 
    battery_charge_pmax = battery_discharge_pmax


    ## Safety net
    # If the available power from the PV is zero at all timesteps, no optimization is needed, since
    # all the consumption is satisfied purchasing energy from the grid. In order to avoid possible
    # problems from the optimization problem in such cases, the optimization is just skipped.

    # Defining a tolerance on the available power to be considered zero (since it is given in kW,
    # a tolerance of 1e-4 is a tenth of a W)
    tol = 1e-4
    if np.all(power_available < 0 + tol):
        # print('Optimization is avoided since there is no excess power from the PV')
        problem_status = 'Opt. unnecessary'
        grid_feed = np.zeros((time_length,))
        grid_purchase = net_load
        battery_charge = np.zeros((time_length,))
        battery_discharge = np.zeros((time_length,))
        battery_energy = battery_energy_min*np.ones((time_length,))

        return problem_status, grid_feed, grid_purchase, battery_charge, battery_discharge, battery_energy



    ### Optimization procedure
    # In case there is excess power from the PV, the optimization procedure is followed

    ## Definition of the problem

    # Initializing the optimization problem using Pulp
    # The problem is set as minimizing the objective function
    opt_problem = LpProblem('Pulp', LpMinimize)
    

    ## Definition of the variables

    # Initalizing the variables
    grid_feed = time_length * [0]  
    grid_feed_state = time_length * [0]   
    grid_purchase = time_length * [0]  
    grid_purchase_state = time_length * [0]  
    battery_charge = time_length * [0]  
    battery_charge_state = time_length * [0]  
    battery_discharge = time_length * [0]  
    battery_discharge_state = time_length * [0]  
    battery_energy = time_length * [0]   

    # Assigning the variables 
    for i in range(time_length):

        # Grid purchase power and state (1|0)        
        if grid_purchase_flag == 1:
            grid_purchase[i] = LpVariable("grid_purchase " + str(i), lowBound=0)
            grid_purchase_state[i] = LpVariable("grid_purchase_state " + str(i), cat=LpBinary)
 
        # Grid feed power and state (1|0)  
        if grid_feed_flag == 1:
            grid_feed[i] = LpVariable("grid_feed " + str(i), lowBound=0)
            grid_feed_state[i] = LpVariable("grid_feed_state " + str(i), cat=LpBinary)
 
        # Battery charge/discharge power and state (1|0)  
        if battery_flag == 1:
            battery_charge[i] = LpVariable("battery_charge " + str(i), lowBound=0)
            battery_charge_state[i] = LpVariable("battery_charge_state " + str(i), cat=LpBinary)
            battery_discharge[i] = LpVariable("battery_discharge " + str(i), lowBound=0)
            battery_discharge_state[i] = LpVariable("battery_discharge_state " + str(i), cat=LpBinary)
            battery_energy[i] = LpVariable("battery_energy " + str(i), lowBound=0)


    ## Constraints (to be set for each time-step)
 
    for i in range(time_length):  

        # Equilibrium at the electric node (in-coming power = exiting power)
        opt_problem += (net_load[i] + grid_feed[i] + battery_charge[i] \
                         - grid_purchase[i] - battery_discharge[i])*dt == 0  #- power_available[i]
        
        # Energy conservation for the battery (and initial SOC = final SOC)
        if (i < time_length - 1):
            opt_problem += (- battery_energy[i + 1] + eta_self_discharge*battery_energy[i]
                            + (battery_charge[i]*eta_charge \
                            - battery_discharge[i]*(1/ eta_discharge))*dt) == 0
        else:
            opt_problem += (- battery_energy[0] + eta_self_discharge*battery_energy[i] \
                            + (battery_charge[i]*eta_charge \
                            - battery_discharge[i]*(1/ eta_discharge))*dt) == 0

        # Constraint on maximum grid power (both for feed and purchase)
        opt_problem += (grid_feed[i] <= grid_feed_state[i] * grid_power_max) 
        opt_problem += (grid_purchase[i] <= grid_purchase_state[i] * grid_power_max)

        # Constraint on feeding/purchasing: they cannot be both active at the same time
        opt_problem += (grid_feed_state[i] + grid_purchase_state[i] >= 0)
        opt_problem += (grid_feed_state[i] + grid_purchase_state[i] <= 1)

        # Constraint on maximum charge and discharge power
        opt_problem += (battery_charge[i] <= battery_charge_state[i] * battery_charge_pmax)
        opt_problem += (battery_discharge[i] <= battery_discharge_state[i] * battery_discharge_pmax)

        # Constraint on charging/discharging: they cannot be both active at the same time
        opt_problem += (battery_charge_state[i] + battery_discharge_state[i] >= 0)
        opt_problem += (battery_charge_state[i] + battery_discharge_state[i] <= 1)

        # Constraint on maximum and minimum SOC
        opt_problem += (battery_energy[i] <= battery_energy_max)
        opt_problem += (battery_energy[i] >= battery_energy_min)

        # Constraint on grid feed: the battery cannot be discharged to sell to the grid
        opt_problem += (grid_feed[i] <= power_available[i]) 

        # Constraint on grid purchase: the battery cannot be charged from the grid
        opt_problem += (battery_charge[i] <= power_available[i]) 
    
    
    ## Setting the objective 
    # The objective is to minimize the interactions with the grid

    opt_problem += lpSum([grid_feed[i] + grid_purchase[i]   for i in range(time_length)])
    
    # The problem is saved in a text file
    with open('opt_problem.txt', 'w') as f:
        print(opt_problem, file=f)   

    
    ## Solution of the problem
    # For each time-step the variables are evaluated in order to reach the objective

    # In some particular cases PULP fails at optimizing the problem and raises an error
    # In order to avoid stopping the procedure due to such errors, a try-except is used
    # If the xception raises, nans are returned

    try:
        opt_problem.solve(PULP_CBC_CMD(msg = 0)) #PULP_CBC_CMD(msg=True)
        
    except:
        optimization_status = 'Opt. did not work'
        grid_feed = np.zeros((time_length,)); grid_feed[:] = np.nan
        grid_purchase = np.zeros((time_length,)); grid_purchase[:] = np.nan
        battery_charge = np.zeros((time_length,)); battery_charge[:] = np.nan
        battery_discharge = np.zeros((time_length,)); battery_discharge[:] = np.nan
        battery_energy = np.zeros((time_length,)); battery_energy[:] = np.nan

        return optimization_status, grid_feed, grid_purchase, battery_charge, battery_discharge, battery_energy      

    # If instead everything goes smooth, the optimization status is printed and the optimized values are returned

    # Optimization status
    optimization_status = LpStatus[opt_problem.status] 

    ## Post-processing
    # The optimized values of the variables are stored in order to be returned

    for i in range(time_length):
        if grid_purchase_flag == 1:
            grid_purchase[i] = value(grid_purchase[i])
            grid_purchase_state[i] = value(grid_purchase_state[i])

        if grid_feed_flag == 1:
            grid_feed[i] = value(grid_feed[i])
            grid_feed_state[i] = value(grid_feed_state[i])

        if battery_flag == 1:
            battery_charge[i] = value(battery_charge[i])
            battery_charge_state[i] = value(battery_charge_state[i])
            battery_discharge[i] = value(battery_discharge[i])
            battery_discharge_state[i] = value(battery_discharge_state[i])
            battery_energy[i] = value(battery_energy[i])


    return optimization_status, \
        np.asarray(grid_feed), np.asarray(grid_purchase), np.asarray(battery_charge), np.asarray(battery_discharge), np.asarray(battery_energy)