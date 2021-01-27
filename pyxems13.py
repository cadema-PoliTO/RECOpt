import numpy as np
import matplotlib.pyplot as plt
import csv
from pulp import *
# import numpy_financial as npf

# funzione che ottimizza l'operazione della batteria, per un giorno di funzionamento

def compute_cons(Pava, Net_load, E_battmax, Nint, time, delta_t): # routine per calcolare E_lgc
        
    # Attivazione Componenti [0/1]  0 Disattivato  1 Attivato
    #E' possibile escludere dall'analisi alcuni componenti del sistema.
    Grid_Purchase=1
    Grid_Sold=1
    Battery=1
     
     ##########################################################################
     # Dati utili
     
    
    # batteria
    file_batt=np.loadtxt("Battery_spec.txt")  
    SOCmin=file_batt[0]
    SOCmax=file_batt[1]
    t_cd_min=file_batt[2]
    eta_bc=file_batt[3]
    eta_bd=file_batt[4]
    eta_sd=file_batt[5]
    
    E_min = SOCmin*E_battmax  
    E_max = SOCmax*E_battmax
    P_bc_max = E_battmax*(1-SOCmin)/t_cd_min  # [kW] charge
    P_bd_max = E_battmax*(1-SOCmin)/t_cd_min  # [kW] discharge
 
       
    # limiti massimi di potenza sulla rete elettrica
    P_smax =50  # [kW] INPUT
    P_pmax = 50  # [kW] INPUT
       
    ###########################################################################
    # Definizione variabili e problema
    # Ottimizzazione tramite pacchetto PULP
    prob = LpProblem('Pulp', LpMinimize)
    
    # variabili
    P_s = len(time) * [0]  
    d_s = len(time) * [0]  
    P_p = len(time) * [0]  
    d_p = len(time) * [0] 
    P_bc = len(time) * [0]  
    P_bd = len(time) * [0]  
    d_bc = len(time) * [0]  
    d_bd = len(time) * [0]  
    E_batt=len(time) * [0] 
       
    for i in range(Nint+1):
        
        if Grid_Purchase == 1:
           d_p[i] = LpVariable("d_p " + str(i), cat=LpBinary)
           P_p[i] = LpVariable("P_p " + str(i), lowBound=0)
 
        if Grid_Sold == 1:
           P_s[i] = LpVariable("P_s " + str(i), lowBound=0)
           d_s[i] = LpVariable("d_s " + str(i), cat=LpBinary)
 
        if Battery == 1:
           P_bc[i] = LpVariable("P_bc " + str(i), lowBound=0)
           d_bc[i] = LpVariable("d_bc " + str(i), cat=LpBinary)
           d_bd[i] = LpVariable("d_bd " + str(i), cat=LpBinary)
           P_bd[i] = LpVariable("P_bd " + str(i), lowBound=0)
           E_batt[i] = LpVariable("E_batt " + str(i), lowBound=0)
     
    ###########################################################################
    # Vincoli
    # Stato finale = Stato iniziale Batterie
    #prob += (SOC_batt[0] - SOC_batt[Nint]== 0)
    #SOC_batt[0]=SOC_ini
    # Bilancio nodo elettrico
              
    for i in range(Nint+1):
        # print(i)    
        prob += (Net_load[i] - P_p[i] + P_s[i] - P_bd[i] + P_bc[i])* delta_t == 0  
        

        if(i<Nint):
           prob +=(-E_batt[i + 1] + eta_sd*E_batt[i] + (P_bc[i] * eta_bc* delta_t - P_bd[i] \
                                                        *(1/ eta_bd)* delta_t))==0
        else:
        #    print('cacaz')
           prob +=(-E_batt[0] + eta_sd*E_batt[i] + (P_bc[i] * eta_bc* delta_t - P_bd[i] \
                                                    *(1/ eta_bd)* delta_t))==0

        prob += (P_s[i] <= d_s[i] * P_smax) # potenza venduta rete
        prob += (P_p[i] <= d_p[i] * P_pmax) # potenza acquistata rete
        prob += (d_p[i] + d_s[i] >= 0)
        prob += (d_p[i] + d_s[i] <= 1)

        prob += (P_bd[i] <= d_bd[i] * P_bd_max)
        prob += (P_bc[i] <= d_bc[i] * P_bc_max)
        prob += (d_bc[i] + d_bd[i] <= 1)
        prob += (d_bc[i] + d_bd[i] >= 0)
        prob += (E_batt[i] <= E_max)
        prob += (E_batt[i] >= E_min)
        prob += (P_s[i] <= Pava[i])  # non puo' scaricare per vendere alla rete
        prob += (P_bc[i] <= Pava[i]) # non puo' caricare dalla rete
        
    ##########################################################################
    # Funzione obiettivo
        
    #economica
    #prob += lpSum([cp[j] * P_p[j] * delta_t - cs[j] * P_s[j] * delta_t for j in range(Nint)])
        
    # minimizzazione scambi con la  rete
    # prob += lpSum([P_s[j]  for j in range(Nint+1)])
    prob += lpSum([P_s[j]+P_p[j]   for j in range(Nint+1)])
    
    

    ###########################################################################
    # Risoluzione
    
    with open('problem.txt', 'w') as f:
        print(prob, file=f)
    status = prob.solve(PULP_CBC_CMD(msg = 0))  #PULP_CBC_CMD(msg=True)
    obj = value(prob.objective)
    # Il programma come impostazione di default appende nelle liste delle variabili precedentemente definite i nomi di queste.
    # Il ciclo for seguente converte i nomi nei corrispettivi valori

    for i in range(Nint+1):
        if Grid_Purchase == 1:
            d_p[i] = value(d_p[i])
            P_p[i] = value(P_p[i])

        if Grid_Sold == 1:
            P_s[i] = value(P_s[i])
            d_s[i] = value(d_s[i])

        if Battery == 1:
            P_bc[i] = value(P_bc[i])
            d_bc[i] = value(d_bc[i])
            d_bd[i] = value(d_bd[i])
            P_bd[i] = value(P_bd[i])
            E_batt[i] = value(E_batt[i])


    # Stato Ottimizzazione
    print("Stato Ottimizzazione:", LpStatus[prob.status] + '\n')
    print('Objective:\n' + str(obj) + '\n')
    
    return P_bc, P_bd, P_s, P_p, E_batt    
  
