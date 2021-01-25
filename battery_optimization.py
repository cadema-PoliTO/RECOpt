import numpy as np
from pulp import *


###############################################################################


# This module contains the method that performs the optimization problem defining
# the operational strategy of the battery


###############################################################################

def battery_optimization(power_available, net_load, battery_capacity, battery_specs, grid_power_max, time_length, time_sim, dt):
    '''
    The method optimizes the operation of the battery during one day.
    Input:

    Output:

    '''

     
    ## Components activation
    # It is possible to exclude some components from the simulation

    grid_purchase_flag = 1
    grid_feed_flag = 1
    battery_flag = 1
     
    
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



    ### Optimization procedure

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
                        - power_available[i] - grid_purchase[i] - battery_discharge[i])*dt == 0  
        
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

    opt_problem += lpSum([grid_feed[j] + grid_purchase[j]   for j in range(time_length)])
    
    
    ## Solution of the problem
    # For each time-step the variables are evaluated in order to reach the objective
    
    # with open('problem.txt', 'w') as f:
    #     print(opt_problem, file=f)
    # status = opt_problem.solve()
    # obj = value(opt_problem.objective)

    opt_problem.solve(PULP_CBC_CMD(msg=0))

    print('\nOptimization status: {}'.format(LpStatus[opt_problem.status]))
    print('\nObjective: {}'.format(value(opt_problem.objective)))

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


    # # Stato Ottimizzazione
    # print("Stato Ottimizzazione:", LpStatus[prob.status] + '\n')
    # print('Objective:\n' + str(obj) + '\n')
    
    return grid_feed, grid_purchase, battery_charge, battery_discharge, battery_energy 