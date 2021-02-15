from PyMouSh import MoulinShape, TimeStamps, Qin_constant, Qin_sinusoidal, Qin_real, calculate_h_S_schoof
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle

#%%
secinday = 24*3600
ZERO_KELVIN = 273.15

head_bf3 = np.load('head_bf3.npy')
temperature_profile=np.array([ZERO_KELVIN, ZERO_KELVIN])


jeme_moulin = pd.read_csv('Field_Data/head_jeme.csv')
jeme_moulin = jeme_moulin.dropna()
h_real = jeme_moulin.head_bed.to_numpy()
t_real = jeme_moulin.soy.to_numpy()

jeme_basin = pd.read_csv('Field_Data/surface_melt_jeme.csv')
jeme_basin = jeme_basin.dropna()
Qin_data = jeme_basin.Qm3s.to_numpy() + 0.1
Qtime_data = jeme_basin.SOY.to_numpy() + 3*3600 #add routing delay -- Jess found 2h not 4. investigate why

time_start = Qtime_data[0]#[int(2*secinday/900)]  
time_end = time_start + 40*secinday #55
timestep = 30*60 #seconds #timestep is longer to reduce the volume of data
time = TimeStamps(time_start,time_end,timestep)

meltwater_input = Qin_real(time, Qin_data, Qtime_data)

#head_L variations
Qin_mean = 0
dQ = 50 
shift1 = 0 * secinday # sin is already shifted from input
variation = Qin_sinusoidal(time,Qin_mean, dQ, shift=shift1)


ice_thickness = 500
initial_head = ice_thickness
channel_length = 500

#head_L variations
# Qin_mean = 300
# dQ = 100 
# shift1 = 0 * secinday # sin is already shifted from input
# head_L= Qin_sinusoidal(time,Qin_mean, dQ, shift=shift1)
# #head_L = 300


moulin = MoulinShape(ice_thickness = ice_thickness,
                      initial_head = initial_head,
                      initial_subglacial_area = np.pi*0.5**2/2, 
                      channel_length = channel_length)

#initialize head
#moulin.head_L=initial_head-50

for idx,t in enumerate(time) :
    #margin moulin
    
    moulin.run1step(time,
                    timestep,
                    meltwater_input[idx],
                    subglacial_baseflow = 0,  
                    head_L = head_bf3[idx+12]-50)#None)#moulin.head_L+variation[idx]) #head_L[idx])
    
picklefile = open('moulin_headL', 'wb')
pickle.dump(moulin, picklefile)
picklefile.close()

#%%


time_year = time/secinday

fig, axs = plt.subplots(3, sharex=False, figsize=(7,7))
fig.tight_layout() 

axs[0].plot(time_year,moulin.dict['meltwater_input_moulin'])
axs[0].set_ylabel('Qin ($m^3/s$)')
axs[0].set_xlim([180,250])

axs[1].plot(t_real/secinday,h_real, color='black', label = 'data')
axs[1].plot(time_year, moulin.dict['head'],label = 'moulin')
axs[1].plot(time_year, moulin.dict['head_L'], label = 'head_L')
axs[1].set_ylabel('Head (m)')
axs[1].set_xlim([180,250])
axs[1].set_ylim([0,500])
axs[1].legend()

axs[2].plot(time_year, moulin.dict['subglacial_radius'],label='moulin')
axs[2].set_ylabel('sub radius (m)')
axs[2].set_xlabel('day of year')
axs[2].set_xlim([180,250])

