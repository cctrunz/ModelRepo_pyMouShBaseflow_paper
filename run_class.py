from PyMouSh import MoulinShape, TimeStamps, Qin_constant, Qin_sinusoidal, Qin_real
import numpy as np

secinday = 24*3600
time_start=0
time_end= 1*secinday
timestep=300
time = TimeStamps(time_start,time_end,timestep)

Qin_mean = 1
dQ = 0.1 
meltwater_input = Qin_sinusoidal(time,Qin_mean, dQ)




moulin = MoulinShape(moulin_radii=1.)

for t in time :
    moulin.run1step(time,meltwater_input)
    print(moulin.head)
    
#%%
import matplotlib.pyplot as plt

plt.figure()
plt.plot(moulin.sim)

#%%

sim = moulin.sim