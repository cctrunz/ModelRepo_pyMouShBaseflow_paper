from PyMouSh import MoulinShape, TimeStamps, Qin_constant, Qin_sinusoidal
import numpy as np
import matplotlib.pyplot as plt
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



# ratio_meltwater = np.ones(len(list_shift))
# ratio_baseflow = np.ones(len(list_shift))

def calc_amplitude(Qin_mean = 0.3, Qin_amplitude = 0.1, bf_mean = 1, bf_amplitude = 0.1):
    head_amplitude = np.ones(len(list_shift))
#    fig, axs = plt.subplots(3, sharex=True, figsize=(8,5))
    
    meltwater_input = Qin_sinusoidal(time, Qin_mean, Qin_amplitude)        
    for idx_shift,shift in enumerate(list_shift):
   
        baseflow = Qin_sinusoidal(time, bf_mean, bf_amplitude, shift=shift*3600) #mean,amplitude
        moulin = MoulinShape(ice_thickness = ice_thickness,
                              initial_head = initial_head,
                              initial_subglacial_area = initial_subglacial_area, 
                              channel_length = channel_length)
        
        for idx,t in enumerate(time) :
        
            moulin.run1step(time,
                            timestep,
                            meltwater_input[idx],
                            overflow = True,
                            subglacial_baseflow = baseflow[idx]) 
            
        idx_nday = 5*24 #days*number of half an hour in a day
        head_portion = moulin.dict['head'][idx_nday:-1]
        head_amplitude[idx_shift] = np.max(head_portion) - np.min(head_portion) #m
        
        # axs[0].plot(moulin.t/secinday, moulin.dict['meltwater_input_moulin'], label='Qin', color='black')
        # axs[0].plot(moulin.t/secinday, moulin.dict['subglacial_baseflow'], label='baseflow')
        # axs[0].set_ylabel('$m^3/s$')
                       
        # axs[1].plot(moulin.t/secinday, moulin.dict['head'])
        # axs[1].set_ylim([0,ice_thickness])
        # axs[1].set_ylabel('Head (m)')
        # axs[1].set_xlabel('days')
        
        # axs[2].plot(moulin.t/secinday, moulin.dict['subglacial_radius'])
        # axs[2].set_ylabel('Sub. Radius (m)')
        # axs[2].set_xlabel('days')
    return head_amplitude
    
#%%

list_shift = np.arange(13) # in hours   

bf_mean      = 0.1 ,0.2, 0.5,   1, #0.5, 0.5, 0.5, 0.5
bf_amplitude = 0.1, 0.1, 0.1, 0.1, #0.2, 0.3, 0.4, 0.5
palette = sns.color_palette('BuPu', n_colors=len(bf_mean))
#sns.palplot(palette)

fig, ax1 = plt.subplots(1, sharex=False, figsize=(5,5))
for idx_bf in np.arange(len(bf_mean)):
    head_amplitude = calc_amplitude(Qin_mean = 0.3, 
                                    Qin_amplitude = 0.1, 
                                    bf_mean =bf_mean[idx_bf], 
                                    bf_amplitude = bf_amplitude[idx_bf])
    ax1.plot(list_shift,(head_amplitude/head_amplitude[0])*100, 
             linestyle='-', linewidth=2, color=palette[idx_bf],
             label='mean = %s $m^3/s$'%bf_mean[idx_bf]) #marker='o',
ax1.set_ylim([0,100])
ax1.set_xlim([0,12])
sns.despine(trim=True)
ax1.set_ylabel('% Amplitude (m)')
ax1.set_xlabel('Lag (h)')
ax1.legend(loc=1, prop={'size': 6})


#%%        
        
list_shift = np.arange(13) # in hours   

bf_mean      = 0.5, 0.5, 0.5, 0.5, 0.5
bf_amplitude = 0.1, 0.2, 0.3, 0.4, 0.5
palette = sns.color_palette('BuPu', n_colors=len(bf_mean))

fig, ax1 = plt.subplots(1, sharex=False, figsize=(5,5))
for idx_bf in np.arange(len(bf_mean)):
    head_amplitude = calc_amplitude(Qin_mean = 0.3, 
                                    Qin_amplitude = 0.1, 
                                    bf_mean =bf_mean[idx_bf], 
                                    bf_amplitude = bf_amplitude[idx_bf])
    ax1.plot(list_shift,(head_amplitude/head_amplitude[0])*100, 
             linestyle='-', linewidth=2, color=palette[idx_bf],
             label='amplitude = %s $m^3/s$'%bf_amplitude[idx_bf]) #marker='o',) #marker='o',
ax1.set_ylim([0,100])
ax1.set_xlim([0,12])
sns.despine(trim=True)
ax1.set_ylabel('% Amplitude (m)')
ax1.set_xlabel('Lag (h)')
ax1.legend(loc=1, prop={'size': 6})


#%%
list_shift = np.arange(13) # in hours   
tmp =  calc_amplitude(Qin_mean = 0.3, 
                                    Qin_amplitude = 0.1, 
                                    bf_mean = 0, 
                                    bf_amplitude = 0)

original_head_amplitude = tmp[0]

bf_mean      = 1, 1,1,1,1
bf_amplitude = 0.1, 0.2, 0.3, 0.4, 0.5

palette = sns.color_palette('BuPu', n_colors=len(bf_mean))

fig, ax1 = plt.subplots(1, sharex=False, figsize=(5,5))
for idx_bf in np.arange(len(bf_mean)):
    head_amplitude = calc_amplitude(Qin_mean = 0.3, 
                                    Qin_amplitude = 0.1, 
                                    bf_mean =bf_mean[idx_bf], 
                                    bf_amplitude = bf_amplitude[idx_bf])
    ax1.plot(list_shift,(head_amplitude/original_head_amplitude)*100, 
             linestyle='-', linewidth=2, color=palette[idx_bf],
             label='amplitude = %s $m^3/s$'%bf_amplitude[idx_bf]) #marker='o',) #marker='o',
#ax1.set_ylim([0,1])
ax1.set_xlim([0,12])
sns.despine(trim=True)
ax1.set_ylabel('% Amplitude (m)')
ax1.set_xlabel('Lag (h)')
ax1.legend(loc=3, prop={'size': 6})










#sns.cubehelix_palette()





# ax3.plot(ratio_baseflow,head_amplitude,marker='o')
# ax3.set_ylabel('Amplitude (m)')
# ax3.set_xlabel('ratio (m)')        



    
    
    
    #plt.savefig('compare_baseflow')





































