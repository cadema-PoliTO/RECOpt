# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 15:35:26 2020

@author: giamm
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import random
import math

import datareader

###############################################################################


# This script is used to plot the results from the simulation


###############################################################################

# The basepath of the file is stored in a variable 
basepath = Path(__file__).parent

# An /Output/Figures folder is created in order to store the graphs as .png files
dirname = 'Output'
subdirname = 'Figures'

try: Path.mkdir(basepath / dirname / subdirname)
except Exception: pass 


## Scale settings

# This is done in order to convert the data (time, powers, energies) from the units of measures
# used in the simulation (min, W, Wh respectively) to the ones used for the post-processing of the results
ts_dict = {'min':1/1,'h':1/60}
ps_dict = {'W':1/1,'kW':1/1e3}
es_dict = {'Wh':1/1,'kWh':1/1e3,'MWh':1/1e6}


## Building seasons and week dictionaries 

# This is done in order to explore all the seasons and, for each season, both 
# types of days (weekday and weekend)
seasons = dict({'winter':(0,'w'),'spring':(1,'ap'),'summer':(2,'s'),'autumn':(3,'ap')})
days = dict({'week-day':(0,'wd'),'weekend-day':(1,'we')})


## Creating a list of different colors

colors = [(230, 25, 75),
        (60, 180, 75),
        (255, 225, 25),
        (0, 130, 200),
        (245, 130, 48),
        (145, 30, 180),
        (70, 240, 240),
        (240, 50, 230),
        (210, 245, 60),
        (250, 190, 212),
        (0, 128, 128),
        (220, 190, 255),
        (170, 110, 40),
        (255, 250, 200),
        (128, 0, 0),
        (170, 255, 195),
        (128, 128, 0),
        (255, 215, 180),
        (0, 0, 128),
        (128, 128, 128)]

# Transforming into rgb triplets
colors_rgb = []
for color in colors:
    color_rgb = []
    for value in color:
        color_rgb.append(float(value)/255) 
    colors_rgb.append(tuple(color_rgb))


## Seasonal load profiles
# A method for plotting the load profiles for both day-types, for each season is created.

def seasonal_load_profiles(time, powers, season, plot_specs, fig_specs, appliances_data, **params):

    ''' The method returns a figure-handle where seasonal load profiles are plotted, for 
        both day types.
    
    Inputs:
        time - 1d array, vector of time 
        powers - 3d array, different types (axis = 0) of time-load profiles to be plotted (axis = 1) for each day-type (axis = 2)
        season - str, containing the name of the season
        plot_specs - dict, for each type of load profile, the type of plot (str, bar plot or plot) and the legend (str)
        fig_specs - dict, containing specific indications such suptitle, etc.
        **params(if not specified, default values are used)
            - 'time_scale': str, 's', 'h'
            - 'power_scale': str,'W', 'kW'
            - 'energy_scale': str, 'Wh', 'kWh', MWh'
            - 'figsize': tup, height and width of the figure
            - 'orientation': str, 'horizontal', 'vertical'
            - 'font_small': float, size of the small fonts (ticks, ...)
            - 'font_medium': float, size of the medium fonts (legend, labels, ...)
            - 'font_large': float, size of the large fonts (titles, ...)
            - all the simulation parameters (n_hh, en_class, etc.)
    
    Outputs:
        fig - figure handle
    '''
   

    ## Input data for the appliances
    # Appliances' attributes, energy consumptions and user's coefficients 

    # # apps is a 2d-array in which, for each appliance (rows) and attribute value is given (columns)
    # apps_ID = appliances_data['apps_ID']

    # # apps_attr is a dictionary in which the name of each attribute (value) is linked to its columns number in apps (key)
    # apps_attr = appliances_data['apps_attr']
       

    ## Parameters
    
    # Default parameters
    def_params = {
    'time_scale': 'h',
    'power_scale': 'kW',
    'energy_scale': 'MWh',
    'figsize': (297/25.4 , 420/25.4),
    'orientation': 'horizontal',
    'font_small': 14,
    'font_medium': 16,
    'font_large': 18,
    }
    
    # The parameters that are not specified when the function is called are set to the default value
    for param in def_params: 
        if param not in params: params[param] = def_params[param]

    ## Updating parameters

    # Scales setup: factors needed to turn the values of time, power and energy in the correct scale
    time_scale = params['time_scale']
    power_scale = params['power_scale']
    # energy_scale = params['energy_scale']

    ts = ts_dict[time_scale]
    ps = ps_dict[power_scale]
    # es = es_dict[energy_scale]

    # Adjusting time, power and energy to the proper scales
    time = time*ts
    powers = powers*ps
    
    ##
    # Figure setup: figure size and orientation, font-sizes 
    figsize = params['figsize']
    orientation = params['orientation']

    if orientation == 'horizontal': figsize = figsize[::-1]

    fontsize_title = params['font_large']
    fontsize_legend = params['font_medium']
    fontsize_labels = params['font_medium']
    # fontsize_text = params['font_medium']
    fontsize_ticks = params['font_small']
    # fontsize_pielabels = params['font_small']

    ##
    # Creating a figure with multiple subplots, with two rows (one for each type of day)
    fig, ax = plt.subplots(2, 1, sharex = False, sharey = False, figsize = figsize)
    
    suptitle = fig_specs['suptitle']
    fig.suptitle(suptitle, fontsize = fontsize_title, fontweight = 'bold')
    fig.subplots_adjust(left = 0.1, bottom = 0.1, right = 0.9, top = 0.85, wspace = None, hspace = 0.3)
   
    ##
    # Evaluating the time-step of the time-vector in order to set the bars' width
    dt = float((time[-1] - time[0])/(np.size(time) - 1))

    # Evaluating the number of profiles passed to the function for each day-type 
    # It is given for ganted that only two day-types are considered
    n_lps = np.size(powers, axis = 0)

    
    ##
    #Running through the day-types (week-day and weekend-day)
    for day in days:

        # Number corresponding to the type of day (0: week-day, 1: week-end -day)
        dd = days[day][0] 

        # Running through the types of load profiles to be plotted for each day-type
        for lp in range(n_lps):

            # Selecting the correct power-data to plot and the plot specifications
            power = powers[lp, :, dd]
            plot_type = plot_specs[lp][0]
            label = plot_specs[lp][1]

            if plot_type == 'plot':
                ax[dd].plot(time + dt/2, power, color = colors_rgb[lp], linestyle = '-', label = label)
 
            elif plot_type == 'bar':
                ax[dd].bar(time, power, color = colors_rgb[lp], width = dt, align = 'edge', label = label)
        
        title = '{}, {}'.format(season.capitalize(), day)
        ax[dd].set_title(title, fontsize = fontsize_title)
        
    ##
    # Making the figure look properly
    for axi in ax.flatten():
        axi.set_xlabel('Time ({})'.format(time_scale), fontsize = fontsize_labels)
        axi.set_ylabel('Power ({})'.format(power_scale), fontsize = fontsize_labels)
        axi.set_xlim([time[0], time[-1]])
        ymin = np.min(powers); ymax = np.max(powers)
        axi.set_ylim([0.9*ymin, 1.1*ymax])
        # Set one tick each hour on the x-axis
        axi.set_xticks(list(time[: : int(60*ts/dt)]))
        axi.tick_params(axis ='both', labelsize = fontsize_ticks)
        axi.tick_params(axis ='x', labelrotation = 0)
        axi.grid()
        axi.legend(loc = 'upper left', fontsize = fontsize_legend, ncol = 2)
    
    ##
    return(fig)


## Seasonal energy consumptions
# A method for plotting the energy consumption by season (total or average) from appliances or classes of appliances 
# is created.

def seasonal_energy(labels_dict, energies, fig_specs, appliances_data, **params):

    ''' The method returns a figure-handle where total energy consumption for appliances or class of appliances
        is plotted, divided by season.
    
    Inputs:
        labels_dict - dict, containing for each label (keys) a unique corresponding index (value)
        energies - 2d array, energy consumption for each appliance/appliance type (axis = 0), by season (axis = 1) 
        season - str, containing the name of the season
        plot_specs - dict, for each type of load profile, the type of plot (str, bar plot or plot) and the legend (str)
        fig_specs - dict, containing specific indications such suptitle, etc.
        **params(if not specified, default values are used)
            - 'time_scale': str, 's', 'h'
            - 'power_scale': str,'W', 'kW'
            - 'energy_scale': str, 'Wh', 'kWh', MWh'
            - 'figsize': tup, height and width of the figure
            - 'orientation': str, 'horizontal', 'vertical'
            - 'font_small': float, size of the small fonts (ticks, ...)
            - 'font_medium': float, size of the medium fonts (legend, labels, ...)
            - 'font_large': float, size of the large fonts (titles, ...)
            - all the simulation parameters (n_hh, en_class, etc.)
    
    Outputs:
        fig - figure handle

    '''


    ## Input data for the appliances
    # Appliances' attributes, energy consumptions and user's coefficients 

    # # apps is a 2d-array in which, for each appliance (rows) and attribute value is given (columns)
    # apps_ID = appliances_data['apps_ID']

    # # apps_attr is a dictionary in which the name of each attribute (value) is linked to its columns number in apps (key)
    # apps_attr = appliances_data['apps_attr']
    

    ## Parameters

    # Default parameters
    def_params = {
    'time_scale': 'h',
    'power_scale': 'kW',
    'energy_scale': 'MWh',
    'figsize': (297/25.4 , 420/25.4),
    'orientation': 'horizontal',
    'font_small': 14,
    'font_medium': 16,
    'font_large': 18,
    }

    # Setting the parameters that are not specified when the function is called to the default value
    for param in def_params: 
        if param not in params: params[param] = def_params[param]

    # Scales setup: factors needed to turn the values of time, power and energy in the correct scale
    # time_scale = params['time_scale']
    # power_scale = params['power_scale']
    energy_scale = params['energy_scale']

    # ts = ts_dict[time_scale]
    # ps = ps_dict[power_scale]
    es = es_dict[energy_scale]

    # Adjusting energies to the proper scales
    energies = energies*es

    # # Making sure that each value for the energy consumption correspond to the correct label
    # id_list = [labels_dict[label][0] for label in labels_dict]
    # energies = energies[id_list,:]

    ##
    # Figure setup: figure size and orientation, font-sizes 
    figsize = params['figsize']
    orientation = params['orientation']

    if orientation == 'horizontal': figsize = figsize[::-1]

    fontsize_title = params['font_large']
    fontsize_legend = params['font_medium']
    fontsize_labels = params['font_medium']
    fontsize_text = params['font_medium']
    fontsize_ticks = params['font_small']
    # fontsize_pielabels = params['font_small']

    ##
    # Creating a figure 
    fig, ax = plt.subplots(figsize=figsize)

    suptitle = fig_specs['suptitle']
    fig.suptitle(suptitle, fontsize = fontsize_title, fontweight = 'bold')
    fig.subplots_adjust(left = 0.1, bottom = 0.2, right = 0.9, top = 0.88, wspace = None, hspace = 0.3)

    ##
    # Labels for the plot, not sorted
    labels_notsort = [label.capitalize().replace('_',' ') for label in labels_dict]

    # Sum of seasonal energy consumptions, for sorting the data 
    total_heights = np.sum(energies, axis = 1)
    ymax = np.max(total_heights)

    # Indices for slicing and sorting the labels and the energies in increasing order
    heights_sortind = np.argsort(total_heights)
    
    # Labels for the plot, sorted
    # labels = labels_notsort[heights_sortind]
    labels = [labels_notsort[ind] for ind in heights_sortind]

    # Initializing the bottoms to zero, in order to make a stack bar plot
    bottoms = np.zeros(len(labels))

    # Initializing the list of seasons, for the legend
    legend = []

    # Initializing the text to add (total seasonal energy consumption),
    text_to_add = 'Energy consumption by season'

    ##
    # Running through the seasons
    for season in seasons:
        
        # Number corresponding to the season (0: winter, 1: summer, 2: spring, 3: autumn)
        ss = seasons[season][0]

        # Energies corresponding to the seasonal consumption, sorted
        heights = energies[heights_sortind, ss]
        
        # Plotting the energy consumption for each season and updating the bottom values
        ax.bar(labels, heights, bottom = bottoms)
        bottoms = bottoms + heights
        
        # Adding the current season to the legend
        legend.append(season.capitalize())

        # Adding the total energy consumption for the current season to the text to be added
        text_to_add = '\n\n'.join((text_to_add, '{0:s}: {1:.2f} {2:s}'.format(season.capitalize(), np.sum(heights), energy_scale)))

    ##
    # Making the figure look properly
    ax.set_ylim([0, 1.1*ymax])
    ax.set_ylabel('Energy consumption ({}/year)'.format(energy_scale), fontsize = fontsize_labels)
    ax.tick_params(axis ='both', labelsize = fontsize_ticks)
    ax.tick_params(axis ='x', labelrotation = 45)
    ax.grid(axis = 'y')
    ax.legend(legend, loc = 'upper left', ncol = len(seasons), fontsize = fontsize_legend)

    # Adding the text with the total energy consumptions by season
    props = dict(boxstyle='square', facecolor = colors_rgb[2], pad = 0.3, alpha = 0.5)                     
    ax.text(0.02, 0.9, text_to_add.rstrip(), fontsize = fontsize_text, ha = 'left', va = 'top', transform = ax.transAxes , bbox = props)

    ##
    return(fig)


## Yearly energy
# A method for plotting the yearly energy consumption (total/average) for appliances or classes of appliances
# is created.

def yearly_energy(labels_dict, energies, fig_specs, appliances_data, **params):

    ''' The method returns a figure-handle where the yearly energy consumption (total/average) 
    for appliances or class of appliances is plotted.
    
    Inputs:
        labels_dict - dict, containing for each label (keys) a unique corresponding index (value)
        energies - 1d array, energy consumption for each appliance/appliance type (axis = 0) for one year
        plot_specs - dict, for each load profile, the type of plot (str, bar plot or plot) and the legend (str)
        season - str, containing the name of the season
        **params(if not specified, default values are used)
            - 'time_scale': str, 's', 'h'
            - 'power_scale': str,'W', 'kW'
            - 'energy_scale': str, 'Wh', 'kWh', MWh'
            - 'figsize': tup, height and width of the figure
            - 'orientation': str, 'horizontal', 'vertical'
            - 'font_small': float, size of the small fonts (ticks, ...)
            - 'font_medium': float, size of the medium fonts (legend, labels, ...)
            - 'font_large': float, size of the large fonts (titles, ...)
    
    Outputs:
        fig - figure handle

    '''
    

    ## Input data for the appliances
    # Appliances' attributes, energy consumptions and user's coefficients 

    # # apps is a 2d-array in which, for each appliance (rows) and attribute value is given (columns)
    # apps_ID = appliances_data['apps_ID']

    # # apps_attr is a dictionary in which the name of each attribute (value) is linked to its columns number in apps (key)
    # apps_attr = appliances_data['apps_attr']

    
    ## Parameters

    # Default parameters
    def_params = {
    'time_scale': 'h',
    'power_scale': 'kW',
    'energy_scale': 'MWh',
    'figsize': (297/25.4 , 420/25.4),
    'orientation': 'horizontal',
    'font_small': 14,
    'font_medium': 16,
    'font_large': 18,
    }

    # Setting the parameters that are not specified when the function is called to the default value
    for param in def_params: 
        if param not in params: params[param] = def_params[param]

    # Scales setup: factors needed to turn the values of time, power and energy in the correct scale
    # time_scale = params['time_scale']
    # power_scale = params['power_scale']
    energy_scale = params['energy_scale']

    # ts = ts_dict[time_scale]
    # ps = ps_dict[power_scale]
    es = es_dict[energy_scale]

    # Adjusting energies to the proper scales
    energies = energies*es

    # # Making sure that each value for the energy consumption correspond to the correct label
    # id_list = [labels_dict[label][0] for label in labels_dict]
    # energies = energies[id_list,:]

    ##
    # Figure setup: figure size and orientatio, font-sizes 
    figsize = params['figsize']
    orientation = params['orientation']

    if orientation == 'horizontal': figsize = figsize[::-1]

    fontsize_title = params['font_large']
    # fontsize_legend = params['font_medium']
    fontsize_labels = params['font_medium']
    fontsize_text = params['font_medium']
    fontsize_ticks = params['font_small']
    # fontsize_pielabels = params['font_small']


    ##
    # Creating a figure 
    fig, ax = plt.subplots(figsize = figsize)

    suptitle = fig_specs['suptitle']
    fig.suptitle(suptitle, fontsize = fontsize_title, fontweight = 'bold')
    fig.subplots_adjust(left = 0.1, bottom = 0.2, right = 0.9, top = 0.88, wspace = None, hspace = 0.3)

    ##
    # Labels for the plot, not sorted
    labels_notsort = [label.capitalize().replace('_',' ') for label in labels_dict]

    # Indices for slicing and sorting the labels and the energies in increasing order
    heights_sortind = np.argsort(energies)

    # Heights for the plot, sorted
    heights = energies[heights_sortind]
    ymax = np.max(heights)
    
    # Labels for the plot, sorted
    # labels = labels_notsort[heights_sortind]
    labels = [labels_notsort[ind] for ind in heights_sortind]
    
    ##
    # Plotting the yearly energy consumption
    ax.bar(labels, heights) 

    # Showing the value of the energy consumption above each bar
    for index, value in enumerate(heights):
        plt.text(index, 1.01*value, '%.f' %(value), ha = 'center', va = 'bottom', rotation = 0, fontsize = fontsize_text)

    ##
    # Evaluating the percentage of the total energy consumption
    total_energy = np.sum(energies)
    energy_perc = energies/total_energy*100
    heights = np.sort(energy_perc)

    # Creating a twin y-axis and plotting the percentage energy consumption
    ax_tw = ax.twinx()
    ax_tw.plot(labels, heights, 'rs')

    if fig_specs['text'] == 'on':
        # Adding the text showing the total yearly energy consumption
        text_toadd = '\n'.join(('Total energy consumption','{} {}/year'.format(total_energy, energy_scale)))
        props = dict(boxstyle = 'square', facecolor = colors_rgb[2], pad = 0.3, alpha = 0.5)                     
        ax.text(0.02, 0.95, text_toadd, fontsize = fontsize_text, ha ='left', va ='top', transform = ax.transAxes, bbox = props)

    ##
    # Making the figure look properly
    ax.set_ylabel('Energy consumption ({}/year)'.format(energy_scale), fontsize=fontsize_labels)
    ax.set_ylim([0, 1.1*ymax])
    ax.tick_params(axis = 'both', labelsize = fontsize_ticks)
    ax.tick_params(axis = 'x', rotation = 70)
    ax.grid(axis = 'y')

    ax_tw.set_ylabel('Energy consumption (%)', fontsize = fontsize_labels)
    ax_tw.yaxis.label.set_color('r')
    ax_tw.set_ylim([0, 1.1*ymax/total_energy*100])

    ax_tw.spines['right'].set_color('r')
    ax_tw.tick_params(axis = 'y', colors = 'r', labelsize = fontsize_ticks)

    ##
    return(fig)


## Seasonal energy pie
# A method for plotting the yearly energy consumption for appliances or classes of appliances (percentage over total)
# as a pie plot is created.

def seasonal_energy_pie(labels_dict, energies, fig_specs, appliances_data, **params):

    ''' The method returns a figure-handle where the percentage over the total energy consumption 
    for appliances or class of appliances is plotted, divided by season.
    
    Inputs:
        labels_dict - dict, containing for each label (keys) a unique corresponding index (value)
        energies - 2d array, energy consumption for each appliance/appliance type (axis = 0), by season (axis = 1) 
        # plot_specs - dict, for each load profile, the type of plot (str, bar plot or plot) and the legend (str)
        **params(if not specified, default values are used)
            - 'time_scale': str, 's', 'h'
            - 'power_scale': str,'W', 'kW'
            - 'energy_scale': str, 'Wh', 'kWh', MWh'
            - 'figsize': tup, height and width of the figure
            - 'orientation': str, 'horizontal', 'vertical'
            - 'font_small': float, size of the small fonts (ticks, ...)
            - 'font_medium': float, size of the medium fonts (legend, labels, ...)
            - 'font_large': float, size of the large fonts (titles, ...)
    
    Outputs:
        fig - figure handle

    '''


    ## Input data for the appliances
    # Appliances' attributes, energy consumptions and user's coefficients 

    # # apps is a 2d-array in which, for each appliance (rows) and attribute value is given (columns)
    # apps_ID = appliances_data['apps_ID']

    # # apps_attr is a dictionary in which the name of each attribute (value) is linked to its columns number in apps (key)
    # apps_attr = appliances_data['apps_attr']
    

    ## Parameters

    # Default parameters
    def_params = {
    'time_scale': 'h',
    'power_scale': 'kW',
    'energy_scale': 'MWh',
    'figsize': (297/25.4 , 420/25.4),
    'orientation': 'horizontal',
    'font_small': 14,
    'font_medium': 16,
    'font_large': 18,
    }
    
    ##
    # Setting the parameters that are not specified when the function is called to the default value
    for param in def_params: 
        if param not in params: params[param] = def_params[param]

    # Scales setup: factors needed to turn the values of time, power and energy in the correct scale
    # time_scale = params['time_scale']
    # power_scale = params['power_scale']
    # energy_scale = params['energy_scale']

    # ts = ts_dict[time_scale]
    # ps = ps_dict[power_scale]
    # es = es_dict[energy_scale]

    # Figure setup: figure size and orientatio, font-sizes 
    figsize = params['figsize']
    orientation = params['orientation']

    if orientation == 'horizontal': figsize = figsize[::-1]

    fontsize_title = params['font_large']
    fontsize_legend = params['font_medium']
    # fontsize_labels = params['font_medium']
    # fontsize_text = params['font_medium']
    # fontsize_ticks = params['font_small']
    fontsize_pielabels = params['font_small']

    ##
    # Creating a new figure, with multiple subplots (one for each season)
    fig, ax = plt.subplots(2, 2, figsize=figsize)

    suptitle = fig_specs['suptitle']
    fig.suptitle(suptitle, fontsize =fontsize_title , fontweight = 'bold')
    fig.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.88, wspace=None, hspace=0.1)

    # Indices that run through the row and columns of ax
    subpl_row = 0
    subpl_col = 0

    ##
    # RUnning through the seasons
    for season in seasons:

        # Number corresponding to the season (0: winter, 1: summer, 2: spring, 3: autumn)
        ss = seasons[season][0]

        # Updating the row and columns indeces 
        if subpl_col > 1: 
            subpl_row = 1
            subpl_col = 0
        
        # Initializing lists where to store the labels, sizes and colors for the pie plot
        labels = []
        sizes = []
        colors = []

        # Selecting the data to plot (only if the seasonal energy consumption is larger than zero)
        for label in labels_dict: 
            labels_index = labels_dict[label][0]

            if energies[labels_index, ss] > 0: 
                labels.append(label.capitalize().replace('_', ' '))
                sizes.append(energies[labels_index, ss])
                colors.append(colors_rgb[labels_index])
        
    
        ##
        # Plotting the data and making the figure look properly
        ax[subpl_row,subpl_col].pie(sizes, autopct='%1.1f%%', pctdistance=0.6, radius=0.8, frame=True, colors = colors , startangle=60, textprops = {'fontsize':fontsize_pielabels})
        ax[subpl_row,subpl_col].legend(labels)
        ax[subpl_row,subpl_col].axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle
        ax[subpl_row,subpl_col].set_title(season.capitalize(), loc='left', pad = 0.5, fontsize = fontsize_legend, fontweight = 'bold')
        ax[subpl_row,subpl_col].set_xticks([])
        ax[subpl_row,subpl_col].set_yticks([])
        
        subpl_col += 1

    ##
    return(fig)
