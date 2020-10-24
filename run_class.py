from PyMouSh import MoulinShape, TimeStamps, Qin_constant, Qin_sinusoidal, Qin_real, calculate_h_S_schoof
import numpy as np

secinday = 24*3600
time_start = 300
time_end = 300#1*secinday
timestep = 300
time = TimeStamps(time_start,time_end,timestep)

Qin_mean = 1
dQ = 0.1 
meltwater_input = Qin_sinusoidal(time,Qin_mean, dQ)




moulin = MoulinShape(moulin_radii=1.)
#%%
for t in time :
    moulin.run1step(time,timestep,meltwater_input)
    print(moulin.head)
    
# #%%
# import matplotlib.pyplot as plt

# sim = moulin.sim

# # plt.figure()
# # plt.plot(sim)

# #%%
# print(sim[0]['head'])


#%%
from scipy.integrate import solve_ivp

dt = moulin.dt
initial_head = 400
initial_subglacial_area = 1
moulin_area = moulin.moulin_area
z = moulin.z
ice_thickness = moulin.ice_thickness
overburden_pressure = moulin.overburden_pressure
channel_length = moulin.channel_length
Qin_compensated = moulin.Qin_compensated
C3 = moulin.C3
overflow = moulin.overflow


sol = solve_ivp(calculate_h_S_schoof,
                        # initial time and end time. We solve for one timestep.
                        [0, dt],
                        # initial head and channel cross-section area. Uses values in previous timestep.
                        [initial_head, initial_subglacial_area],
                        # args() only works with scipy>=1.4. if below, then error message: missing 5 arguments
                        args=(moulin_area,z, ice_thickness, overburden_pressure,channel_length ,Qin_compensated, C3, overflow),
                        method='LSODA'  # solver method
                        # atol = 1e-6, #tolerance. can make it faster
                        # rtol = 1e-3,
                        # max_step = 10 #change if resolution is not enough
                        )
# (m) moulin water level
head = sol.y[0][-1]
# (m) Channel cross-section
subglacial_area = sol.y[1][-1]

