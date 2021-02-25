from pyMouSh.pyMouSh import MoulinShape, TimeStamps, Qin_constant, Qin_sinusoidal
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
#import pandas as pd
import pickle
import seaborn as sns
sns.set()
sns.set_theme(style="ticks", font_scale=1.25)

secinday = 24*3600
#ZERO_KELVIN = 273.15

dz = 1
z_elevations = None
#moulin_radii = 0.3
#temperature_profile = np.array([ZERO_KELVIN, ZERO_KELVIN])
ice_thickness = 500
initial_head = 500
initial_subglacial_area = np.pi*0.8**2/2      
regional_surface_slope = 0
channel_length = 25e3


time_start = 0
time_end = time_start + 10*secinday #55
timestep = 30*60 #seconds #timestep is longer to reduce the volume of data
time = TimeStamps(time_start,time_end,timestep)

idx_nday = 5*24 #days*number of half an hour in a day

from numpy.random import default_rng
expected_head_amplitude = 60 # based of visual measurements of the head amplitude in JEME data
original_head_amplitude = 532.44 

n_run = 100
# create random numbers 
rng1 = default_rng()
rgn2 = default_rng() 
Qin_mean = 1 * rng1.random(n_run)
Qin_amplitude = rgn2.random(n_run)


#initialize variables for loop:
head_amplitude_exp = np.zeros(n_run)
head_amplitude_ori = np.zeros(n_run)
head = []
meltwater_input = []
idx_plot = []




for idx_run in np.arange(n_run):   
    melt_in = Qin_sinusoidal(time, Qin_mean[idx_run], Qin_amplitude[idx_run])         
    moulin = MoulinShape(ice_thickness = ice_thickness,
                          initial_head = initial_head,
                          initial_subglacial_area = initial_subglacial_area, 
                          channel_length = channel_length)
    if Qin_mean[idx_run]>=Qin_amplitude[idx_run]:    #prevent head
        idx_plot.append(idx_run)
        meltwater_input.append(melt_in)
        
        for idx,t in enumerate(time) :
            moulin.run1step(time,
                            timestep,
                            melt_in[idx],
                            overflow = True,
                            subglacial_baseflow = 0) 
            
        head.append(moulin.dict['head'])
        head_portion = moulin.dict['head'][idx_nday:-1]
        #calculate head amplitude
        head_amplitude_exp[idx_run] = np.max(head_portion) - np.min(head_portion) #m
        head_amplitude_ori[idx_run] = np.max(head_portion) - np.min(head_portion) #m
    



#%%Plots
purple = (0.5326874279123414, 0.2502883506343714, 0.612641291810842)
palette = sns.color_palette('BuPu', n_colors=n_run) #RdBu PuOr
        
fig, ax = plt.subplots(1, sharex=False, figsize=(5,5))     
ax.plot(Qin_amplitude/Qin_mean,(head_amplitude_exp/original_head_amplitude)*100, 
                     linestyle='',marker='o', color=purple)
ax.plot([0,0.35],[100,100],linestyle='--', color='grey')
ax.plot([0.35,0.35],[0,100],linestyle='--', color='grey')
# ax.plot([0.35],[0],marker='v', color='grey')
# ax.plot([0.02],[100],marker='<', color='grey')

#ax.set_ylim([0,800])
ax.set_xlim([0,1])
ax.set_ylabel('% Initial head amplitude (m)')
ax.set_xlabel('meltwater input variability (-)')#('bf_amplitude/bf_mean')
sns.despine(offset=[0,-12.4],trim=True)

plt.savefig('sensitivity_to_meltwater_variability.pdf')

#%%
fig, axs = plt.subplots(2, sharex=True, figsize=(5,5))   
for idx in np.arange(len(idx_plot)):
          
    axs[0].plot(time/secinday, meltwater_input[idx], label='Qin')#, color=palette[idx])
    axs[0].set_ylabel('$m^3/s$')
                   
    axs[1].plot(time/secinday, head[idx])#, color=palette[idx])
    axs[1].set_ylim([0,ice_thickness])
    axs[1].set_ylabel('Head (m)')
    axs[1].set_xlabel('days')
sns.despine(trim=True)



    
