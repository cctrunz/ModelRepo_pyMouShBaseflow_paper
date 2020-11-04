from PyMouSh import MoulinShape, TimeStamps, Qin_sinusoidal, Qin_real
import numpy as np
import matplotlib.pyplot as plt
from celluloid import Camera
import pandas as pd

secinday = 24*3600
ZERO_KELVIN = 273.15


temperature_profile=np.array([ZERO_KELVIN, ZERO_KELVIN])
ice_thickness = 500
regional_surface_slope = 0
channel_length = 15000

jeme_basin = pd.read_csv('Field_Data/surface_melt_jeme.csv')
jeme_basin = jeme_basin.dropna()
Qin_data = jeme_basin.Qm3s.to_numpy() +0.1
Qtime_data = jeme_basin.SOY.to_numpy()

jeme_moulin = pd.read_csv('Field_Data/head_jeme.csv')
jeme_moulin = jeme_moulin.dropna()
h_real = jeme_moulin.head_bed.to_numpy()
t_real = jeme_moulin.soy.to_numpy()

time_start = Qtime_data[0]
time_end = Qtime_data[0] + 10*secinday# Qtime_data[-1]#time_start + 5*secinday #
timestep = 300 #seconds
time = TimeStamps(time_start,time_end,timestep)

meltwater_input = Qin_real(time, Qin_data, Qtime_data)
# Qin_mean = 1
# dQ = 0.1 
# meltwater_input = Qin_sinusoidal(time,Qin_mean, dQ)
Q=pd.DataFrame(meltwater_input)
Q.rolling(200)
subglacial_baseflow = Q.to_numpy()


sim = MoulinShape(  
                    channel_length = channel_length,
                    temperature_profile = temperature_profile,                   
                    ice_thickness = ice_thickness,
                    regional_surface_slope = regional_surface_slope)

for t in time :
    #main subglacial channel
    sim.run1step(time,
                    timestep,
                    meltwater_input,
                    subglacial_baseflow = subglacial_baseflow)

fig = plt.figure(figsize=(7,6),dpi=100) 
sim.plot_AGU_3(fig,-2,t_real,h_real)   
# PLOT MODEL OUTPUT   
# plt.figure()
# idx = int(moulin_baseflow.idx/2)
# moulin_baseflow.plot_AGU(idx,t_real,h_real,spine_head_min=200)
#%%
# fig, ax = plt.subplots()
# moulin_baseflow.plot_head(ax)

#fig, ax = plt.subplots()
#sim.plot_Qin(ax)
#sim.plot_radius(ax,z_position = 'heq',bottom_axis=True)
fig = plt.figure(figsize=(7,6),dpi=100)
#fig, ax, = plt.subplots()
camera = Camera(fig)
for idx in np.arange(0,sim.idx,12):
    #sim.plot_radius(ax)
    #sim.plot_head(ax,idx_max=idx)
    sim.plot_AGU_3(fig,idx,t_real,h_real)
    camera.snap()

animation = camera.animate(interval = 100)#, repeat = True,repeat_delay = 500)
animation.save('testagu3.mp4')#,dpi=200)




# fig = plt.figure()
# camera = Camera(fig)
# for i in range(10):
#     plt.plot([i] * 10)
#     camera.snap()
# animation = camera.animate()
# animation.save('celluloid_minimal.gif', writer = 'imagemagick')