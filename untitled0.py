from PyMouSh import MoulinShape, TimeStamps, Qin_sinusoidal
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import pandas as pd
import os

secinday = 24*3600
ZERO_KELVIN = 273.15

channel_length = 30000
ice_thickness = 1000

temperature_profile =np.array([ZERO_KELVIN, ZERO_KELVIN])
regional_surface_slope = 0
initial_subglacial_area = 1#(np.pi*0.2**2)/2)

#time
start = 0
end = 5
timestep = 300
time = np.arange(start*secinday,end*secinday,timestep)
time_plot = np.arange(start*secinday,end*secinday,timestep*12)


#idealized Qin for test purposes
Qin_mean = 3
dQ = 0.1 
meltwater_input = Qin_sinusoidal(time,Qin_mean, dQ)

moulin = MoulinShape( moulin_radii = 4.,
                     initial_head = 800,
                      channel_length = channel_length,
                      ice_thickness = ice_thickness)
                     
for idx,t in enumerate(time):
    moulin.run1step(t,timestep, meltwater_input[idx],
                    creep=True,
                    elastic_deformation=False,
                    melt_below_head=False,
                    open_channel_melt=False,
                    potential_drop=False,
                    ice_motion=False,
                    refreezing=False)
 
#%%


for idx,t in enumerate(time_plot):
    plt.figure(figsize=(2,5))
    if idx == 0:
        plt.plot(moulin.listdict[idx]['delta_creep_major']*1000,moulin.z,color='firebrick')        
        
    if idx == 1:
        plt.plot(moulin.listdict[idx-1]['delta_creep_major']*1000,moulin.z,color='r')
        plt.plot(moulin.listdict[idx]['delta_creep_major']*1000,moulin.z,color='firebrick')
        
    if idx == 2:
        plt.plot(moulin.listdict[idx-2]['delta_creep_major']*1000,moulin.z,color='tomato')
        plt.plot(moulin.listdict[idx-1]['delta_creep_major']*1000,moulin.z,color='r')
        plt.plot(moulin.listdict[idx]['delta_creep_major']*1000,moulin.z,color='firebrick')
        
    if idx == 3:
        plt.plot(moulin.listdict[idx-3]['delta_creep_major']*1000,moulin.z,color='lightsalmon')
        plt.plot(moulin.listdict[idx-2]['delta_creep_major']*1000,moulin.z,color='tomato')
        plt.plot(moulin.listdict[idx-1]['delta_creep_major']*1000,color='r')
        plt.plot(moulin.listdict[idx]['delta_creep_major']*1000,moulin.z,color='firebrick')
        
    if idx >3:
        plt.plot(moulin.listdict[idx-4]['delta_creep_major']*1000,moulin.z,color='peachpuff')
        plt.plot(moulin.listdict[idx-3]['delta_creep_major']*1000,moulin.z,color='lightsalmon')
        plt.plot(moulin.listdict[idx-2]['delta_creep_major']*1000,moulin.z,color='tomato')
        plt.plot(moulin.listdict[idx-1]['delta_creep_major']*1000,moulin.z,color='r')
        plt.plot(moulin.listdict[idx]['delta_creep_major']*1000,moulin.z,color='firebrick')
        
    #plt.xlim(-10,0)
    

    