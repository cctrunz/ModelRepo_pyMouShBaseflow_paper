#TEST OC

from PyMouSh import MoulinShape, TimeStamps, Qin_constant, Qin_sinusoidal
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#%%
secinday = 24*3600
ZERO_KELVIN = 273.15

time_start = 0
time_end = 5*secinday
timestep = 300 #seconds
time = TimeStamps(time_start,time_end,timestep)

Qin_mean = 3
dQ = 0.1 
meltwater_input = Qin_sinusoidal(time,Qin_mean, dQ)


#paramters                          
friction_factor_OC = 1
friction_factor_TM = 0.5
friction_factor_SUB = 0.1


moulin = MoulinShape(dz=1,                    
                    friction_factor_OC = friction_factor_OC,
                    friction_factor_TM = friction_factor_TM,
                    friction_factor_SUB = friction_factor_SUB)

for t in time :
    #main subglacial channel
    moulin.run1step(time,
                    timestep,
                    meltwater_input,
                    subglacial_baseflow = 0, 
                    head_L = None,
                    open_channel_melt=True,
                    potential_drop=False)

#idx = -2
fig, ax = plt.subplots()
moulin.plot_moulin(ax)