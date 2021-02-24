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

tmp = pd.read_csv('Field_Data/temperature_foxx1.csv')
temperature_profile = tmp.temperature.to_numpy()
#temperature_profile =np.array([ZERO_KELVIN, ZERO_KELVIN])
regional_surface_slope = 0.1
initial_subglacial_area = 1#(np.pi*0.2**2)/2)

#time
time = 300
timestep = 300


meltwater_input = 3
moulin_radii = 4.
initial_head = 400#np.arange(400,1000,10)


moulin = MoulinShape( moulin_radii = moulin_radii,
                      initial_head = initial_head,
                      channel_length = channel_length,
                      ice_thickness = ice_thickness)
                     

moulin.run1step(time,timestep, meltwater_input,
                creep=True,
                elastic_deformation=True,
                melt_below_head=True,
                open_channel_melt=True,
                potential_drop=False,
                ice_motion=True,
                refreezing=False)
 
#%%


plt.plot(moulin.listdict[0]['delta_creep_major']*1000,moulin.z,color='firebrick')        
plt.plot(moulin.listdict[0]['delta_elastic_major']*1000,moulin.z,color='orange')   
plt.plot(moulin.listdict[0]['delta_melt_below_head']*1000,moulin.z,color='dodgerblue')            
plt.plot(moulin.listdict[0]['delta_melt_above_head_open_channel']*1000,moulin.z,color='dodgerblue') 
plt.plot(moulin.listdict[0]['delta_ice_motion']*1000,moulin.z,color='mediumvioletred')       

    #plt.xlim(-10,0)

#%%    

    