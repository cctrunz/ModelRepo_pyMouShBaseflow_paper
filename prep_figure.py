import numpy as np

#from typing import Union, Tuple#, ArrayLike
import matplotlib.pyplot as plt
from collections import defaultdict

from pyMouSh import MoulinShape, TimeStamps, Qin_constant, Qin_sinusoidal, Qin_real, calculate_h_S_schoof
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle


secinday = 24*3600
ZERO_KELVIN = 273.15
#paramters

#temperature profile based of foxx data -- Lauren
tmp = pd.read_csv('Field_Data/temperature_foxx1.csv')
temperature_profile = tmp.temperature.to_numpy()#np.array([ZERO_KELVIN, ZERO_KELVIN])

moulin_radii = 0.3
initial_subglacial_area = 1 
channel_length = 25e3
ice_thickness = 500
initial_head = ice_thickness

#Import meltwater input calculated with meltmodel for JEME
jeme_basin = pd.read_csv('Field_Data/surface_melt_jeme.csv')
jeme_basin = jeme_basin.dropna()
Qin_data = jeme_basin.Qm3s.to_numpy() +0.1
Qtime_data = jeme_basin.SOY.to_numpy() + 3*3600 #add routing delay -- Jess found 2h not 4. investigate why

#Import head measurement for JEME
jeme_moulin = pd.read_csv('Field_Data/head_jeme.csv')
jeme_moulin = jeme_moulin.dropna()
h_real = jeme_moulin.head_bed.to_numpy()
t_real = jeme_moulin.soy.to_numpy()

#set timeserie length and timestep for interpolation
time_start = Qtime_data[int(2*secinday/900)]  
time_end = time_start + 50*secinday
timestep = 30*60#300 #seconds to reduce number of saved data
time = TimeStamps(time_start,time_end,timestep)

#calculate meltwater input (interpolate)
meltwater_input = Qin_real(time, Qin_data, Qtime_data)

#calculate baseflow
bf_mean = 1
bf_amp = 0.1 
shift0h = 0.42 * secinday
shift12h = 0 * secinday # sin is already shifted from input
baseflow_shift0 = Qin_sinusoidal(time,bf_mean, bf_amp, shift=shift0h)
baseflow_shift12 = Qin_sinusoidal(time,bf_mean, bf_amp, shift=shift12h)




def find_nearest(array, value):
    """Finds the nearest value in an array and outputs a index.
    This function was found in 
    https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array

    Parameters:
    -----------
    array: array to be looked into
    value: single value to look for into the array

    Output:
    -------
    index of the closest value in the array
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


def plot_radius(sim,axis,idx_max=-2, z_position='heq', bottom_axis=True):
     """ Plot the moulin radii timeserie for a certain position in the moulin. 
     
     Parameters:
     -----------
         axis (string) : ax .. what is that, how to describe it??
         z_position (int, string) : altitude in the moulin where the moulin radius will be display.
                                     'heq': at equilibrium head
         bottom_axis (optional, bolean) : True --> display the time axis
                                          False --> don't display time axis
     
     Example:
     --------
         fig, ax = plt.subplots()
         sim.plot_head(ax)
         """

     Mr_z = defaultdict(list)

     if z_position == 'heq':
         for idx in sim.dict['all_idx']:
             idx_position = find_nearest(sim.z,sim.dict['head'][idx])  
             Mr_z['major'].append(sim.listdict[idx]['moulin_radius_major'][idx_position])
             Mr_z['minor'].append(sim.listdict[idx]['moulin_radius_minor'][idx_position])         
     else:
         idx_position = find_nearest(sim.z,z_position)
         for idx in sim.dict['all_idx']: 
             Mr_z['major'].append(sim.listdict[idx]['moulin_radius_major'][idx_position])
             Mr_z['minor'].append(sim.listdict[idx]['moulin_radius_minor'][idx_position])
         
     axis.plot(sim.time_day[0,idx_max],Mr_z['major'][0,idx_max],'-',color='black') 
     axis.plot(sim.time_day[0,idx_max],Mr_z['minor'][0,idx_max],'-',color='grey') 
     
     min_time = int(min(sim.time_day))
     max_time = int(max(sim.time_day))
     
     axis.set_ylabel('(m)',color='black')
     axis.set_xlim([min_time,max_time])
     axis.set_ylim([min(Mr_z['major']),max(Mr_z['major'])]) 
     axis.yaxis.tick_left()
     axis.yaxis.set_label_position("left")
     axis.tick_params(axis='y', labelcolor='black')
     axis.spines['top'].set_visible(False)
     
     axis.spines['right'].set_visible(False)
     axis.spines['left'].set_color('black')
     if bottom_axis == False:
         axis.spines['bottom'].set_visible(False)
         axis.axes.xaxis.set_visible(False)
     #axis.spines['left'].set_bounds()
     #axis.set_yticks(np.arange()) 
     #axis.set_yticklabels(np.arange())

     
def plot_moulin(sim, axis, idx, 
                 left_lim = -5,
                 left_bound = -4,
                 right_lim = 5,
                 right_bound = 4,
                 ground_depth = -50,
                 x_tick_spacing = 2,
                 axis_side = 'left'
                 ):
     """ Plot moulin shape and head in the moulin for a certain timestep 
     Parameters:
     -----------
         axis
         idx (optional, int): index run to display
     """
     
     Mx_upstream = sim.listdict[idx]['moulin_wall_position_upstream']
     Mx_downstream = sim.listdict[idx]['moulin_wall_position_downstream']      
     
     axis.plot(Mx_upstream[0:-10],sim.z[0:-10],color='black') #plot major axis on the left
     axis.plot(Mx_downstream[0:-10],sim.z[0:-10],color='black')  #plot minor axis on the right

     axis.axhspan(0, sim.dict['head'][idx], facecolor ='lightblue', alpha = 1,zorder=1)
     axis.axhspan(ground_depth, 0, facecolor ='peru', alpha = 1,zorder=1)
     axis.fill_betweenx(sim.z,left_lim,Mx_upstream,color='aliceblue',zorder=2)
     axis.fill_betweenx(sim.z,Mx_downstream,right_lim,color='aliceblue',zorder=2)      

     axis.set_ylabel('Elevation (m)')    
     axis.set_xlabel('(m)')
     axis.set_xlim([left_lim,right_lim]) 
     axis.set_ylim([ground_depth,sim.ice_thickness]) 
     axis.spines['top'].set_visible(False)
     axis.spines['bottom'].set_position(('zero')) 
     axis.spines['bottom'].set_bounds(right_bound,left_bound)
     axis.set_xticks(np.arange(left_bound,right_bound+1,x_tick_spacing))
     axis.set_xticklabels(abs(np.arange(left_bound,right_bound+1,x_tick_spacing)) )
     axis.set_yticks(np.arange(0,sim.ice_thickness+1,100)) 
     axis.set_yticklabels(np.arange(0,sim.ice_thickness+1,100))  
     
     if axis_side == 'right':
         axis.yaxis.tick_right()
         axis.yaxis.set_label_position('right')        
         axis.spines['left'].set_visible(False)
         axis.spines['right'].set_bounds(0,sim.ice_thickness)  
                   
     if axis_side == 'left':
         axis.yaxis.tick_left()
         axis.yaxis.set_label_position('left')        
         axis.spines['right'].set_visible(False)
         axis.spines['left'].set_bounds(0,sim.ice_thickness)   

def plot_head(sim,axis,
               idx_min=0,
               idx_max=-1,
               spine_head_min=0,
               bottom_axis=True,
               color='black',
               axis_side = 'left',
               ground_depth=-50):
     """ Plot head timeserie 
     
     Parameters:
     -----------
         axis (string) : ax .. what is that, how to describe it??
         spine_head_min (optional, int) : min value to display 
                 on the y axis on the graph
                 default value is zero (display all the way to the bed)
         bottom_axis (optional, bolean) : True --> display the time axis
                                          False --> don't display time axis
     
     Example:
     --------
         fig, ax = plt.subplots()
         MoulinShape.plot_head(ax)
         """       

     axis.plot(sim.time_day[idx_min:idx_max],sim.dict['head'][idx_min:idx_max],'-',color=color)        
     axis.set_ylabel('$m$',color=color,rotation='horizontal')           
     axis.set_xlim([sim.time_day[idx_min],sim.time_day[idx_max]])
     axis.set_ylim([ground_depth,sim.ice_thickness]) 
     
     axis.tick_params(axis='y', labelcolor=color)
     axis.spines['top'].set_visible(False)
     
     axis.set_yticks(np.arange(spine_head_min,sim.ice_thickness+1,100)) 
     axis.set_yticklabels(np.arange(spine_head_min,sim.ice_thickness+1,100))
     
     if axis_side == 'right':
         axis.yaxis.tick_right()
         axis.yaxis.set_label_position('right')        
         axis.spines['left'].set_visible(False)
         axis.spines['right'].set_color(color)
         axis.tick_params(axis='y', colors=color)
         axis.spines['right'].set_bounds(spine_head_min,sim.ice_thickness)
         
     if axis_side == 'left':
         axis.yaxis.tick_left()
         axis.yaxis.set_label_position('left')               
         axis.spines['right'].set_visible(False)
         axis.spines['left'].set_color(color)
         axis.tick_params(axis='y', colors=color)
         axis.spines['left'].set_bounds(spine_head_min,sim.ice_thickness)
     
     if bottom_axis == False:
         axis.spines['bottom'].set_visible(False)
         axis.axes.xaxis.set_visible(False)

def plot_subglacial_radius(sim,axis,
                            idx_min=0,
                            idx_max=-1,
                            bottom_axis=True,
                            color='black',
                            axis_side = 'left'):
     '''Subglacial channel'''
     axis.plot(sim.time_day[idx_min:idx_max],sim.dict['subglacial_radius'][idx_min:idx_max],'-',color=color)  
     axis.set_ylabel('$m$',color=color,rotation='horizontal')   
     axis.set_xlim([sim.time_day[idx_min],sim.time_day[idx_max]])
     axis.set_ylim([min(sim.dict['subglacial_radius']),max(sim.dict['subglacial_radius'])])
     axis.tick_params(axis='y', labelcolor=color)
     axis.spines['top'].set_visible(False)
     
     if axis_side == 'right':
         axis.yaxis.tick_right()
         axis.yaxis.set_label_position('right')        
         axis.spines['left'].set_visible(False)
         axis.spines['right'].set_color(color)
         axis.tick_params(axis='y', colors=color)
         
     if axis_side == 'left':
         axis.yaxis.tick_left()
         axis.yaxis.set_label_position('left')        
         axis.spines['right'].set_visible(False)
         axis.spines['left'].set_color(color)
         axis.tick_params(axis='y', colors=color)
         
     if bottom_axis == False:
         axis.spines['bottom'].set_visible(False)
         axis.axes.xaxis.set_visible(False)
     
def plot_subglacial_channel(sim,axis,
                             idx_min=0,
                             idx_max=-1,
                             bottom_axis=True,
                             color='black',
                             axis_side = 'left'   
                             ):
     '''Subglacial channel'''
     axis.plot(sim.time_day[idx_min:idx_max],sim.dict['subglacial_cross_section_area'][idx_min:idx_max],'-',color=color)  
     axis.set_ylabel('$m^2$',color=color,rotation='horizontal')   
     #axis.set_xlim([sim.time_day[idx_min],sim.time_day[idx_max]])
     axis.set_ylim([min(sim.dict['subglacial_cross_section_area']),max(sim.dict['subglacial_cross_section_area'])])
     axis.tick_params(axis='y', labelcolor=color)
     axis.spines['top'].set_visible(False)
     
     if axis_side == 'right':
         axis.yaxis.tick_right()
         axis.yaxis.set_label_position('right')        
         axis.spines['left'].set_visible(False)
         axis.spines['right'].set_color(color)
         axis.tick_params(axis='y', colors=color)
         
     if axis_side == 'left':
         axis.yaxis.tick_left()
         axis.yaxis.set_label_position('left')        
         axis.spines['right'].set_visible(False)
         axis.spines['left'].set_color(color)
         axis.tick_params(axis='y', colors=color)
         
     if bottom_axis == False:
         axis.spines['bottom'].set_visible(False)
         axis.axes.xaxis.set_visible(False)
     
def plot_Qin(sim,axis,
              idx_min=0,
              idx_max=-1,
              bottom_axis=True,
              color='black',
              axis_side = 'left'
              ):        
     '''Qin'''
     axis.plot(sim.time_day[idx_min:idx_max],sim.dict['meltwater_input_moulin'][idx_min:idx_max],'-',color=color)#,label='Qin') 
     axis.set_ylabel('$m^3/s$',color=color,rotation='horizontal') 
     axis.set_xlim([sim.time_day[idx_min],sim.time_day[idx_max]])
     axis.set_ylim([min(sim.dict['meltwater_input_moulin']),max(sim.dict['meltwater_input_moulin'])])
     axis.spines['top'].set_visible(False)
     axis.tick_params(axis='y', labelcolor=color)       
     
     if axis_side == 'right':
         axis.yaxis.tick_right()
         axis.yaxis.set_label_position('right')        
         axis.spines['left'].set_visible(False)
         axis.spines['right'].set_color(color)
         axis.tick_params(axis='y', colors=color)
         
     if axis_side == 'left':
         axis.yaxis.tick_left()
         axis.yaxis.set_label_position('left')        
         axis.spines['right'].set_visible(False)
         axis.spines['left'].set_color(color)
         axis.tick_params(axis='y', colors=color)
     if bottom_axis == False:
         axis.spines['bottom'].set_visible(False)
         axis.axes.xaxis.set_visible(False)   
         
     #ax4.set_xticks(np.round(np.linspace(min_time,20,max_time)))
     #ax4.set_xticklabels(np.round(np.linspace(min_time,20,max_time)))

def plot_baseflow(sim,axis,
                   idx_min=0,
                   idx_max=-1,
                   bottom_axis=True,
                   color='black',
                   axis_side = 'left'):        
     '''Baseflow'''
     axis.plot(sim.time_day[idx_min:idx_max],sim.dict['subglacial_baseflow'][idx_min:idx_max],'-',color=color)#,label='Qin') 
     axis.set_ylabel('$m^3/s$',color=color,rotation='horizontal')   
     axis.set_xlim([sim.time_day[idx_min],sim.time_day[idx_max]])
     #axis.set_ylim([min(sim.dict['subglacial_baseflow']),max(sim.dict['subglacial_baseflow'])])
     axis.spines['top'].set_visible(False)
     axis.tick_params(axis='y', labelcolor=color)        
     
     if axis_side == 'right':
         axis.yaxis.tick_right()
         axis.yaxis.set_label_position('right')        
         axis.spines['left'].set_visible(False)
         axis.spines['right'].set_color(color)
         axis.tick_params(axis='y', colors=color)
         
     if axis_side == 'left':
         axis.yaxis.tick_left()
         axis.yaxis.set_label_position('left')        
         axis.spines['right'].set_visible(False)
         axis.spines['left'].set_color(color)
         axis.tick_params(axis='y', colors=color)
         
     if bottom_axis == False:
         axis.spines['bottom'].set_visible(False)
         axis.axes.xaxis.set_visible(False) 

     
def plot_Qout(sim,axis,
               idx_min=0,
               idx_max=-1,
               bottom_axis=True,
               color='black',
               axis_side = 'left'):        
     '''Qout'''
     axis.plot(sim.time_day[idx_min:idx_max],sim.dict['meltwater_output_subglacial'][idx_min:idx_max],'-',color=color)#,label='Qin') 
     axis.set_ylabel('$m^3/s$',color=color,rotation='horizontal')           
     axis.set_xlim([sim.time_day[idx_min],sim.time_day[idx_max]])
     axis.set_ylim([min(sim.dict['meltwater_output_subglacial']),max(sim.dict['meltwater_output_subglacial'])])
     axis.tick_params(axis='y', labelcolor=color)
     axis.spines['top'].set_visible(False)
     
     if axis_side == 'right':
         axis.yaxis.tick_right()
         axis.yaxis.set_label_position('right')        
         axis.spines['left'].set_visible(False)
         axis.spines['right'].set_color(color)
         axis.tick_params(axis='y', colors=color)
         
     if axis_side == 'left':
         axis.yaxis.tick_left()
         axis.yaxis.set_label_position('left')        
         axis.spines['right'].set_visible(False)
         axis.spines['left'].set_color(color)
         axis.tick_params(axis='y', colors=color)
         
     if bottom_axis == False:
         axis.spines['bottom'].set_visible(False)
         axis.axes.xaxis.set_visible(False)   
         
         
         
def plot_data(moulin_sim,t_real,h_real,idx=-1):
    
    fig = plt.figure(figsize=(13,5))   
    grid = plt.GridSpec(4,3)#, wspace=-0.7)
    
    ax1b = fig.add_subplot(grid[0, 0:2])#Qin
    ax1a = ax1b.twinx() #baseflow 
    ax2 = fig.add_subplot(grid[1:4, 2])#moulin
    ax3 = fig.add_subplot(grid[1:4, 0:2])#hw
    ax4 = fig.add_subplot(grid[3, 0:2])#SCs  

    #ax1a.yaxis.set_label_coords(-0.09, 0.5)
    ax1b.yaxis.set_label_coords(-0.065, 0.5)
    ax3.yaxis.set_label_coords(-0.08, 0.7)
    ax4.yaxis.set_label_coords(-0.08, 0.5)
    
    ground_depth=-60
    spine_head_min=200
    lw=1.5
    
    #baseflow
    moulin_sim.plot_baseflow(ax1a,
                       color='seagreen',
                       bottom_axis=False,
                       axis_side = 'right') 
    #ax1a.set_xlim(t_lim) 

    #Meltwater
    moulin_sim.plot_Qin(ax1b,
                   bottom_axis=False,
                   axis_side = 'left',
                   color='grey') 
    #ax1b.set_xlim(t_lim) 
    #ax1b.set_ylim(Q_lim)


    #Moulin
    moulin_sim.plot_moulin(ax2,idx,
                     left_lim = -5,
                     left_bound = -4,
                     right_lim = 5,
                     right_bound = 4,
                     x_tick_spacing = 2,
                     ground_depth=ground_depth,
                     axis_side = 'right',)

    #Head
    ax3.plot(t_real/3600/24,h_real,'-',color='black') 
    moulin_sim.plot_head(ax3,
                   color='steelblue',
                   spine_head_min = spine_head_min,
                   bottom_axis = False,
                   axis_side = 'left',
                   ground_depth = ground_depth) 

    #ax3.set_xlim(t_lim)


    #Subglacial radius
    moulin_sim.plot_subglacial_radius(ax4,
                               color='orangered',
                               bottom_axis=True,
                               axis_side = 'left')         
    #ax4.set_xlim(t_lim)
    #ax4.tick_params(axis='x',labelsize=8)
    #ax4.set_ylim(SC_lim)
    #ax4.spines['bottom'].set_bounds(t_lim[0]+0.1,t_lim[1])
    ax4.set_xlabel('Day of year 2017')


    #Legend  
    l1a = ax1a.legend(['Subglacial baseflow'],loc="upper left", frameon=True, bbox_to_anchor=(0.6, 1.1))
    #loc=1, prop={'size': 12},bbox_to_anchor=(1.1, 1.6), bbox_transform=axs[1].transAxes)
    
    for line, text in zip(l1a.get_lines(), l1a.get_texts()):
        text.set_color(line.get_color())  
    l1b = ax1b.legend(['Meltwater input'],loc="upper left", bbox_to_anchor=(0, 1.1) )
    for line, text in zip(l1b.get_lines(), l1b.get_texts()):
        text.set_color(line.get_color())  
    l3 = ax3.legend(['Head measured','Head simulated'],loc="upper left", bbox_to_anchor=(0, 1) )    
    for line, text in zip(l3.get_lines(), l3.get_texts()):
        text.set_color(line.get_color())
    l4 = ax4.legend(['Subglacial radius'],loc="upper left", bbox_to_anchor=(0, 1.1) )
    for line, text in zip(l4.get_lines(), l4.get_texts()):
        text.set_color(line.get_color())
    ax4.patch.set_alpha(0)

    #thicker axis
    for axis in ['top','bottom','left','right']:
        ax1a.spines[axis].set_linewidth(lw)
    ax1a.tick_params(width=lw)
    
    for axis in ['top','bottom','left','right']:
        ax1b.spines[axis].set_linewidth(lw)
    ax1b.tick_params(width=lw)

    for axis in ['top','bottom','left','right']:
        ax2.spines[axis].set_linewidth(lw)
    ax2.tick_params(width=lw)

    for axis in ['top','bottom','left','right']:
        ax3.spines[axis].set_linewidth(lw)
    ax3.tick_params(width=lw)

    for axis in ['top','bottom','left','right']:
        ax4.spines[axis].set_linewidth(lw)
    ax4.tick_params(width=lw)
        

#%%

bf0_fix = MoulinShape(moulin_radii = moulin_radii,
                      ice_thickness = ice_thickness,
                      initial_head = initial_head,
                      initial_subglacial_area = initial_subglacial_area,
                      channel_length = channel_length,
                      temperature_profile = temperature_profile)

for idx,t in enumerate(time) :
    bf0_fix.run1step(t,
                    timestep,
                    meltwater_input[idx],
                    subglacial_baseflow = 0)
plot_data(bf0_fix,t_real,h_real) 