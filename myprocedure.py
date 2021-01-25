import numpy as np
#from numpy import *
import matplotlib.pyplot as plt
from pulp import *
import csv
# import numpy_financial as npf
from tictoc import tic, toc
from scipy.interpolate import interp1d
from aggregated_households_load_profiler import aggregate_households_load_profiler as agr_hlp


def battery_optimization(power_available, net_load, battery_capacity, battery_specs, time_length, time_sim, dt):
        
    # Attivazione Componenti [0/1]  0 Disattivato  1 Attivato
    #E' possibile escludere dall'analisi alcuni componenti del sistema.
    grid_purchase_flag = 1
    grid_feed_flag = 1
    battery_flag = 1
     
     ##########################################################################
     # Dati utili

    
     
    
    # batteria
    
    SOCmin = battery_specs['SOC_min']
    SOCmax = battery_specs['SOC_max']
    t_cd_min = battery_specs['t_cd_min']
    eta_charge = battery_specs['eta_charge']
    eta_discharge = battery_specs['eta_discharge']
    eta_self_discharge = battery_specs['eta_self_discharge']
    
    battery_energy_min = SOCmin*battery_capacity 
    battery_discharge_pmax = battery_capacity*(1-SOCmin)/t_cd_min  # [kW] discharge

    print(t_cd_min)
    print(battery_discharge_pmax)

    battery_energy_max = SOCmax*battery_capacity 
    battery_charge_pmax = battery_discharge_pmax
      
    # limiti massimi di potenza sulla rete elettrica
    grid_power_max = 50  # [kW] INPUT
    grid_power_min = 50  # [kW] INPUT
       
    ###########################################################################
    # Definizione variabili e problema
    # Ottimizzazione tramite pacchetto PULP
    opt_problem = LpProblem('Pulp', LpMinimize)
    
    # variabili
    grid_feed = time_length * [0]  
    grid_feed_state = time_length * [0]   
    grid_purchase = time_length * [0]  
    grid_purchase_state = time_length * [0]  
    battery_charge = time_length * [0]  
    battery_charge_state = time_length * [0]  
    battery_discharge = time_length * [0]  
    battery_discharge_state = time_length * [0]  
    battery_energy = time_length * [0]   
       
    for i in range(time_length):
        
        if grid_purchase_flag == 1:
            grid_purchase[i] = LpVariable("grid_purchase " + str(i), lowBound=0)
            grid_purchase_state[i] = LpVariable("grid_purchase_state " + str(i), cat=LpBinary)
 
        if grid_feed_flag == 1:
            grid_feed[i] = LpVariable("grid_feed " + str(i), lowBound=0)
            grid_feed_state[i] = LpVariable("grid_feed_state " + str(i), cat=LpBinary)
 
        if battery_flag == 1:
            battery_charge[i] = LpVariable("battery_charge " + str(i), lowBound=0)
            battery_charge_state[i] = LpVariable("battery_charge_state " + str(i), cat=LpBinary)
            battery_discharge[i] = LpVariable("battery_discharge " + str(i), lowBound=0)
            battery_discharge_state[i] = LpVariable("battery_discharge_state " + str(i), cat=LpBinary)
            battery_energy[i] = LpVariable("battery_energy " + str(i), lowBound=0)
     
    ###########################################################################
    # Vincoli
    # Stato finale = Stato iniziale Batterie
    #prob += (SOC_batt[0] - SOC_batt[Nint]== 0)
    #SOC_batt[0]=SOC_ini
    # Bilancio nodo elettrico
              
    for i in range(time_length):
        print(i)    
        opt_problem += (net_load[i] + grid_feed[i] + battery_charge[i] \
                        - power_available[i] - grid_purchase[i] - battery_discharge[i])*dt == 0  
        

        if (i < time_length - 1):
            opt_problem += (- battery_energy[i + 1] + eta_self_discharge*battery_energy[i]
                            + (battery_charge[i]*eta_charge \
                            - battery_discharge[i]*(1/ eta_discharge))*dt) == 0
        else:
            print('cielostorto')
            opt_problem += (- battery_energy[0] + eta_self_discharge*battery_energy[i] \
                            + (battery_charge[i]*eta_charge \
                            - battery_discharge[i]*(1/ eta_discharge))*dt) == 0

        opt_problem += (grid_feed[i] <= grid_feed_state[i] * grid_power_max) # potenza venduta rete
        opt_problem += (grid_purchase[i] <= grid_purchase_state[i] * grid_power_max) # potenza acquistata rete
        opt_problem += (grid_feed_state[i] + grid_purchase_state[i] >= 0)
        opt_problem += (grid_feed_state[i] + grid_purchase_state[i] <= 1)

        opt_problem += (battery_charge[i] <= battery_charge_state[i] * battery_charge_pmax)
        opt_problem += (battery_discharge[i] <= battery_discharge_state[i] * battery_discharge_pmax)
        opt_problem += (battery_charge_state[i] + battery_discharge_state[i] >= 0)
        opt_problem += (battery_charge_state[i] + battery_discharge_state[i] <= 1)
        opt_problem += (battery_energy[i] <= battery_energy_max)
        opt_problem += (battery_energy[i] >= battery_energy_min)
        opt_problem += (grid_feed[i] <= power_available[i])  # non puo' scaricare per vendere alla rete
        opt_problem += (battery_charge[i] <= power_available[i]) # non puo' caricare dalla rete
        
    ##########################################################################
    # Funzione obiettivo
        
    #economica
    #prob += lpSum([cp[j] * P_p[j] * delta_t - cs[j] * P_s[j] * delta_t for j in range(Nint)])
        
    # minimizzazione scambi con la  rete
    # opt_problem += lpSum([grid_feed[j]  for j in range(time_length)])
    opt_problem += lpSum([grid_feed[j] + grid_purchase[j]   for j in range(time_length)])
    
    

    ###########################################################################
    # Risoluzione
    
    with open('problem.txt', 'w') as f:
        print(opt_problem, file=f)
    status = opt_problem.solve()
    obj = value(opt_problem.objective)
    # Il programma come impostazione di default appende nelle liste delle variabili precedentemente definite i nomi di queste.
    # Il ciclo for seguente converte i nomi nei corrispettivi valori

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








pv_size = 10                # Potenza di picco PV (kW))
battery_size = 5          # Capacità batteria (kWh)

pv_production_unit = np.loadtxt("PPV.txt")                  # Produzione fotovoltaica (per 1 kWp)
pv_production = pv_production_unit*pv_size                  # Produzione fotovoltaica scalata alla taglia
    
       
ue_seasons = agr_hlp()                   # Carico elettrico delle utenze nella configurazione (kW)

dt = 1                                                      # Time-step (h)
time = 24                                                   # Total time of simulation for each day (h)
time_sim = np.arange(0, time, dt)                           # Vector of time for the simulation (h)
time_length = np.size(time_sim)

seasons_n = 4
months_n = 12 
days_n = 2

seasons = np.arange(0, months_n, int(months_n/seasons_n))    
months = np.arange(0, months_n)
days = np.arange(0, days_n)

print(seasons)
print(months)

ue_consumption = np.zeros((months_n, time_length, days_n))

for day in days:

    # f = interp1d(seasons, ue_seasons[:, :, day], kind = 'linear', axis = 0, fill_value = 'extrapolate')
    # ue_consumption[:, :, day] = f(months)

    for timestep in time_sim:
        ue_consumption[:, timestep, day] = np.interp(months, seasons, ue_seasons[:, timestep, day], period = months_n)

    
    
plt.figure()
plt.plot(time_sim, ue_seasons[0, :, 0], linestyle = '--', label = 'winter')
plt.plot(time_sim, ue_seasons[1, :, 0], linestyle = '--', label = 'spring')
# plt.plot(time_sim, ue_seasons[2, :, 0], linestyle = '--', label = 'summer')
# plt.plot(time_sim, ue_seasons[3, :, 0], linestyle = '--', label = 'autumn')

for month in months[:6]:

    if month not in seasons:
        plt.plot(time_sim, ue_consumption[month, :, 0], label = month)

plt.legend()
plt.show()

# ## Battery specifications
# battery_capacity = battery_size
# battery_specifications_file = np.loadtxt("Battery_spec.txt")  
# battery_specs = {
#     'SOC_min': battery_specifications_file[0],
#     'SOC_max': battery_specifications_file[1],
#     't_cd_min': battery_specifications_file[2],
#     'eta_charge': battery_specifications_file[3],
#     'eta_discharge': battery_specifications_file[4],
#     'eta_self_discharge': battery_specifications_file[5],
#     }

 

# ## Powers to be considered
# # Initializing the vectors
# pv = np.zeros((time_length, months_n))                      # Power from PV that is not directly self-consumed (kW)
# pv_available = np.zeros((time_length, months_n))            # Power from PV tht is not directy or indirectly self-consumed (kW)
# ue_net_load = np.zeros((time_length, months_n))             # Net load for the users part of the configuration (kW)
# ue_sf_net_load = np.zeros((time_length, months_n))          # Net load for the consumption unit directly connectd to the production plant (kW)
# battery_charge= np.zeros((time_length, months_n))           # Power that charges the batter (kW)
# battery_discharge = np.zeros((time_length, months_n))       # Power discharged from the battery (kW)
# grid_feed = np.zeros((time_length, months_n))               # Power fed into the grid (kW)
# grid_purchase = np.zeros((time_length, months_n))           # Power purchased from the grid (kW)
# battery_energy = np.zeros((time_length, months_n))          # Energy stored in the battery (kWh)

# # Evaluating the actual values

# pv = pv_production - ue_sf_consumption
# ue_sf_net_load[pv < 0] = -pv[pv < 0]
# pv[pv < 0] = 0


# pv_available = pv - ue_consumption
# ue_net_load[pv_available < 0] = -pv_available[pv_available < 0]
# pv_available[pv_available < 0] = 0

# net_load = ue_net_load + ue_sf_net_load

# for jj in range(months_n):
#     for ii in range(time_length):

#         if pv_available[ii,jj] > 0 and ue_net_load[ii, jj] > 0: print('There\'s a problem')


# ## Showing some data

# month = 1
# plt.subplot()
# plt.plot(time_sim, ue_sf_consumption[:, month], 'b--', label = 'ue_sf')
# plt.plot(time_sim, ue_consumption[:, month], 'b', label = 'ue')
# plt.plot(time_sim, pv_production[:, month], 'r--', label = 'pv_prod')
# plt.plot(time_sim, pv[:, month], 'r', label = 'pv')
# plt.plot(time_sim, pv_available[:, month], 'g', label = 'pv_ava')
# plt.plot(time_sim, net_load[:, month], 'g--', label = 'net_load')

# plt.legend()
# plt.show()


# ## Solving
# # Mixed Integer Linear Problem (MILP) Optmiziation for the battery control

# # import pyxems13
# if month >= 0 & month < 12:

#     j = month
        
#     grid_feed[:,j], grid_purchase[:,j], \
#     battery_charge[:,j], battery_discharge[:,j], \
#     battery_energy[:,j] = battery_optimization(pv_available[:,j], net_load[:,j], battery_capacity, battery_specs, time_length, time_sim, dt)
    
#     # battery_charge[:,j], battery_discharge[:,j],\
#     # grid_feed[:,j], grid_purchase[:,j], \
#     # battery_energy[:,j] = pyxems13.compute_cons(pv_available[:,j], net_load[:,j], battery_capacity, time_length - 1, time_sim + 1, dt)

# # # Calcolo potenza autoconsumata istantaneamente (condivisa dalla comunità)
# # P_lgc=np.zeros( (Nint+1,Nmonth) )

# # print(battery_energy[:,0])


# for i in range(time_length):

#     if ((net_load[i,month] + grid_feed[i,month] + battery_charge[i,month] \
#                         - pv_available[i,month] - grid_purchase[i,month] - battery_discharge[i,month])*dt) > 1e-6: print('oh oh')


# print(battery_discharge[:, month])
# print(battery_charge[:, month])
# print(grid_feed[:, month])
# print(grid_purchase[:, month])

# fig = plt.figure()
# ax = fig.add_subplot()
# # ax.plot(time_sim, battery_discharge[:, month], 'b--', label = 'P_b_d')
# # ax.plot(time_sim, battery_charge[:, month], 'b', label = 'P_b_c')
# ax.plot(time_sim, grid_purchase[:, month], 'r--', label = 'P_g_p')

# ax.plot(time_sim, ue_consumption[:, month] + ue_sf_consumption[:, month], 'b', label = 'ue')
# # ax.plot(time_sim, grid_feed[:, month], 'r', label = 'P_g_f')
# # plt.plot(time_sim[:, month], 'g', label = 'pv_ava')
# # plt.plot(net_load[:, month], 'g--', label = 'net_load')
# plt.legend()

# ax_twin = ax.twinx()
# ax_twin.plot(time_sim, battery_energy[:, month], 'y', label = 'E_b')
# plt.legend()

# plt.show()



# energy_shared = np.minimum((pv + battery_discharge), (ue_consumption + battery_charge))


# print(np.sum(grid_feed[:, month]))
# print(np.sum(energy_shared[:, month]))
# print(np.sum(energy_shared[:, month])/np.sum(pv[:, month] + battery_discharge[:, month])*100)