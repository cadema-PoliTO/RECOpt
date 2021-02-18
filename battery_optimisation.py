import numpy as np
from pulp import *


###############################################################################


# This module contains the method that performs the optimization problem defining
# the operational strategy of the battery


###############################################################################

def battery_optimization(pv_production, consumption, time_dict, technologies_dict):
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

    # PV size (kW)
    pv_size = technologies_dict['pv_size']

    # Grid maximum power (kW)
    grid_feed_max = technologies_dict['pv_size']
    grid_purchase_max = technologies_dict['grid_power_max']

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
    battery_discharge_max = battery_capacity*(SOCmax-SOCmin)/t_cd_min 
    battery_charge_max = battery_discharge_max

    # # Evaluating the excess energy, if the production tries to fulfill the whole demand (needed to bound the grid_feed power)
    # pv_available = pv_production - consumption
    # pv_available[pv_available < 0] = 0


    # ## Easy solution
    #
    # # If the available power from the PV is zero at all timesteps, no optimization is needed, since
    # # all the consumption is satisfied purchasing energy from the grid. In such cases, the optimization 
    # # can just be skipped.
    #
    # # Defining a tolerance on the available power to be considered zero (since it is given in kW,
    # # a tolerance of 1e-4 is a tenth of a W)
    # tol = 1e-4
    # if np.all(power_available < 0 + tol):
    #     # print('Optimization is avoided since there is no excess power from the PV')
    #     problem_status = 'Opt. unnecessary'
    #     grid_feed = np.zeros((time_length,))
    #     grid_purchase = net_load
    #     battery_charge = np.zeros((time_length,))
    #     battery_discharge = np.zeros((time_length,))
    #     battery_energy = battery_energy_min*np.ones((time_length,))

    #     return problem_status, grid_feed, grid_purchase, battery_charge, battery_discharge, battery_energy



    ### Optimization procedure
    # In case there is excess power from the PV, the optimization procedure is followed
    # Please select the objective to pursue, between:
    # 'MINGRI': MINimize GRid Interactions (grid feed of excess energy and grid purchase of deficit energy)
    # 'MAXSHE': MAXimize SHared Energy (hourly minimum between all energy fed into the grid, PV + battery,
    #           and all energy taken from the grid, consumption + battery)
    opt_objectives = ['MAXSHE', 'MINGRI']

    opt_objective_num = 0
    opt_objective = opt_objectives[opt_objective_num]


    ## Definition of the problem

    # Initializing the optimization problem using Pulp
    # The problem is set as minimizing the objective function

    if opt_objective == 'MINGRI': opt_problem = LpProblem('Pulp', LpMinimize)
    elif opt_objective == 'MAXSHE': opt_problem = LpProblem('Pulp', LpMaximize)
    

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
    shared_power = time_length * [0] 
    
    # y is a binary variable that is used to linearize the definition of the shared power
    # The function min(x1, x2) is indeed not linear therefore a proper implementation is needed
    # y is used to assess whether x1 is larger than x2 or the other way around
    y = time_length * [0]

    # M is a big-parameter that is used in the linearization of the function min(x1, x2)
    # In this case, x1 and x2 are, respectively, the hourly sum between the pv production and the battery discharge,
    # and the hourly sum between the consumption and the battery charge. M is a big number that must be larger than both x1 and x2
    # The battery_charge/discharge are not known in advance but their maximum value is, therefore M is evaluated as follows
    # M = 100*max(np.max(pv_production) + battery_discharge_pmax, \
    #             np.max(consumption) + battery_charge_pmax)
    M = 100*max(pv_size + battery_discharge_max, \
                grid_purchase_max + battery_charge_max)

    # Assigning the variables 
    for i in range(time_length):

        # Grid purchase power and state (1|0)        
        grid_purchase[i] = LpVariable("grid_purchase " + str(i), lowBound = 0)
        grid_purchase_state[i] = LpVariable("grid_purchase_state " + str(i), cat = LpBinary)
 
        # Grid feed power and state (1|0)  
        grid_feed[i] = LpVariable("grid_feed " + str(i), lowBound = 0)
        grid_feed_state[i] = LpVariable("grid_feed_state " + str(i), cat = LpBinary)
 
        # Battery charge/discharge power and state (1|0) and battery energy
        battery_charge[i] = LpVariable("battery_charge " + str(i), lowBound = 0)
        battery_charge_state[i] = LpVariable("battery_charge_state " + str(i), cat = LpBinary)
        battery_discharge[i] = LpVariable("battery_discharge " + str(i), lowBound = 0)
        battery_discharge_state[i] = LpVariable("battery_discharge_state " + str(i), cat = LpBinary)
        battery_energy[i] = LpVariable("battery_energy " + str(i), lowBound = 0)
        
        # Shared power and linearization variable y (1|0)
        shared_power[i] = LpVariable("shared_power " + str(i), lowBound = 0)
        y[i] = LpVariable("auxiliary " + str(i), cat = LpBinary)


    ## Constraints (to be set for each time-step)
 
    for i in range(time_length):  

        # Equilibrium at the electric node (in-coming power = exiting power)
        opt_problem += (consumption[i] + grid_feed[i] + battery_charge[i] \
                         - pv_production[i] - grid_purchase[i] - battery_discharge[i])*dt == 0  #- power_available[i]
        
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
        opt_problem += (grid_feed[i] <= grid_feed_state[i] * grid_feed_max) 
        opt_problem += (grid_purchase[i] <= grid_purchase_state[i] * grid_purchase_max)

        # Constraint on feeding/purchasing: they cannot be both active at the same time
        opt_problem += (grid_feed_state[i] + grid_purchase_state[i] >= 0)
        opt_problem += (grid_feed_state[i] + grid_purchase_state[i] <= 1)

        # Constraint on maximum charge and discharge power
        opt_problem += (battery_charge[i] <= battery_charge_state[i] * battery_charge_max)
        opt_problem += (battery_discharge[i] <= battery_discharge_state[i] * battery_discharge_max)

        # Constraint on charging/discharging: they cannot be both active at the same time
        opt_problem += (battery_charge_state[i] + battery_discharge_state[i] >= 0)
        opt_problem += (battery_charge_state[i] + battery_discharge_state[i] <= 1)

        # Constraint on maximum and minimum SOC
        opt_problem += (battery_energy[i] <= battery_energy_max)
        opt_problem += (battery_energy[i] >= battery_energy_min)

        # Constraint on grid feed: the battery cannot be discharged to sell to the grid
        opt_problem += (grid_feed[i] <= pv_production[i]) 

        # Constraint on grid purchase: the battery cannot be charged from the grid
        opt_problem += (battery_charge[i] <= pv_production[i]) 

        # Linearization of shared_power[i] = min(pv_production[i] + battery_discharge[i], consumption[i] + battery_charge[i])

        # Constraint on the shared energy, that must be smaller than both pv_production + battery_discharge
        # and consumption + battery_charge (1/2)
        opt_problem += (shared_power[i] <= pv_production[i] + battery_discharge[i])
        opt_problem += (shared_power[i] <= consumption[i] + battery_charge[i])

        # Definition of y that is 1 when pv_production + battery_discharge <= consumption[i] + battery_charge[i], 0 otherwise
        # The definition of y is introduced as a constraint
        opt_problem += ((consumption[i] + battery_charge[i]) - (pv_production[i] + battery_discharge[i]) <= M*y[i])
        opt_problem += ((pv_production[i] + battery_discharge[i]) - (consumption[i] + battery_charge[i]) <= M*(1 - y[i]))

        # Constraint on the shared energy, that must be not only smaller than both (...) but also equal to the minimum value
        # when y == 1, shared_power = pv_production[i] + battery_discharge[i] since it is both larger-equal for this constraint
        # and smaller-equal for the previous one. When y == 0, the other way around
        opt_problem += (shared_power[i] >= pv_production[i] + battery_discharge[i] - M*(1 - y[i]))
        opt_problem += (shared_power[i] >= consumption[i] + battery_charge[i]  - M*y[i])
  
    
    ## Setting the objective 

    # Objective of minimizing the interactions with the grid
    if opt_objective == 'MINGRI': opt_problem += lpSum([grid_feed[i] + grid_purchase[i]   for i in range(time_length)])

    # Objective of maximizing the shared energy
    elif opt_objective == 'MAXSHE': opt_problem += lpSum([shared_power[i] for i in range(time_length)])
    

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
        
        # Grid purchase
        grid_purchase[i] = value(grid_purchase[i])
        grid_purchase_state[i] = value(grid_purchase_state[i])

        # Grid feed
        grid_feed[i] = value(grid_feed[i])
        grid_feed_state[i] = value(grid_feed_state[i])

        # Battery charge/ discharge/ energy
        battery_charge[i] = value(battery_charge[i])
        battery_charge_state[i] = value(battery_charge_state[i])
        battery_discharge[i] = value(battery_discharge[i])
        battery_discharge_state[i] = value(battery_discharge_state[i])
        battery_energy[i] = value(battery_energy[i])

        # Shared power
        shared_power[i] = value(shared_power[i])


    return optimization_status, \
        np.asarray(shared_power), np.asarray(grid_feed), np.asarray(grid_purchase), \
        np.asarray(battery_charge), np.asarray(battery_discharge), np.asarray(battery_energy)