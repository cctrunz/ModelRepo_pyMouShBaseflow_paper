from PyMouSh import MoulinShape, TimeStamps, Qin_constant, Qin_sinusoidal
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
initial_subglacial_area = np.pi*0.88**2/2      
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

# n_run = 10
# # create random numbers 
# rng1 = default_rng()
# rgn2 = default_rng() 
# bf_mean      = 3 * rng1.random(n_run) #3*np.ones(n_run)
# bf_amplitude = 1*rgn2.random(n_run)
bf_mean = 1
bf_amplitude = 0.1

#initialize variables for loop:
lags = np.arange(13)
head_amplitude= np.zeros(len(lags))

head = []
baseflow = []

Qin_mean = 0.3
Qin_amplitude = 0.1
meltwater_input = Qin_sinusoidal(time, Qin_mean, Qin_amplitude)


bf_mean      = 0.5
bf_amplitude = 0.1
palette_lag = sns.color_palette('PuOr', n_colors=len(lags))

#%%

fig, axs = plt.subplots(2, sharex=True, figsize=(5,5))   
axs[0].plot(time/secinday, meltwater_input, label='Qin', color='black')

for idx_lag,lag in enumerate(lags):
    bf = Qin_sinusoidal(time, bf_mean, bf_amplitude,shift=lag*3600) #mean,amplitude
    moulin = MoulinShape(ice_thickness = ice_thickness,
                          initial_head = initial_head,
                          initial_subglacial_area = initial_subglacial_area, 
                          channel_length = channel_length)
        
    for idx,t in enumerate(time) :
        moulin.run1step(time,
                        timestep,
                        meltwater_input[idx],
                        overflow = True,
                        subglacial_baseflow = bf[idx])         
    # head.append(moulin.dict['head'])
    head_portion = moulin.dict['head'][idx_nday:-1]
    #calculate head amplitude
    head_amplitude[idx_lag] = np.max(head_portion) - np.min(head_portion) #m
    
    
    axs[0].plot(time/secinday, bf, label='baseflow', color=palette_lag[idx_lag])
    axs[0].set_ylabel('$m^3/s$')
                   
    axs[1].plot(time/secinday, moulin.dict['head'], color=palette_lag[idx_lag])
    axs[1].set_ylim([0,ice_thickness])
    axs[1].set_ylabel('Head (m)')
    axs[1].set_xlabel('days')
    sns.despine(trim=True)
    
plt.savefig('lag_for_all.pdf')    
    

#%%
fig, ax = plt.subplots(1, sharex=False, figsize=(5,5))   

bf_mean      = 0.1 ,0.2, 0.5,   1, #0.5, 0.5, 0.5, 0.5
bf_amplitude = 0.1, 0.1, 0.1, 0.1, #0.2, 0.3, 0.4, 0.5
palette_lag = sns.color_palette('BuPu', n_colors=len(lags))

for idx_run in np.arange(len(bf_mean)):
    for idx_lag,lag in enumerate(lags):
        bf = Qin_sinusoidal(time, bf_mean[idx_run], bf_amplitude[idx_run],shift=lag*3600) #mean,amplitude
        moulin = MoulinShape(ice_thickness = ice_thickness,
                              initial_head = initial_head,
                              initial_subglacial_area = initial_subglacial_area, 
                              channel_length = channel_length)
            
        for idx,t in enumerate(time) :
            moulin.run1step(time,
                            timestep,
                            meltwater_input[idx],
                            overflow = True,
                            subglacial_baseflow = bf[idx])         
        # head.append(moulin.dict['head'])
        head_portion = moulin.dict['head'][idx_nday:-1]
        #calculate head amplitude
        head_amplitude[idx_lag] = np.max(head_portion) - np.min(head_portion) #m
        

#Plots
    ax.plot(lags,(head_amplitude/head_amplitude[0])*100, 
             linestyle='-', linewidth=2,
             label='mean = %s $m^3/s$'%bf_mean[idx_run]) #marker='o', color=palette[idx_run],
    
    
            
ax.set_ylim([0,100])
ax.set_xlim([0,12])
sns.despine(trim=True)
ax.set_ylabel('% total amplitude (m)')
ax.set_xlabel('Lag (h)')
ax.legend(loc=3, prop={'size': 10})

plt.savefig('sensitivity_to_lag.pdf')


    
