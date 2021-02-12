## Compare Moush (0D), fixed cylinder (0D), and fixed cylinder (1D)
import numpy as np
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt


#for the 1D simulation. This requires to have Matt's modules installed
from conduits1D_landlab_matt_02112021 import run1Dsim,plot_3panels_singlestep, plot_3panels_movie, plot_2panel_overtime_multiposition #for the 1D simulation. This requires to have Matt's modules installed

from PyMouSh import MoulinShape, TimeStamps, Qin_sinusoidal, Qin_real

secinday = 24*3600
ZERO_KELVIN = 273.15

#Moulin parameters for MouSh
dz = 1
z_elevations = None
temperature_profile=np.array([ZERO_KELVIN, ZERO_KELVIN])
friction_factor_OC = 1
friction_factor_TM = 0.5
fraction_pd_melting = 0.2
creep_enhancement_factor = 5
sigma = (0, 0)  # -50e3 #(Units??) compressive #50e3#-50e3
tau_xy = 0  # -50e3#100e3 #(Units??) shear opening

#Glacier parameters for JEME
moulin_radii = 0.3 #m (based out of the highest radius when running moush for JEME)
ice_thickness = 500 #m
initial_head = 500 #m
initial_subglacial_area = 1 #m
channel_length = 25e3 #m 
friction_factor_SUB = 0.1 # its the same value in pressurized_flow in landlab
regional_surface_slope = 0

#Import meltwater input for JEME from melt model (from Jessica Mejia)
jeme_basin = pd.read_csv('Field_Data/surface_melt_jeme.csv')
jeme_basin = jeme_basin.dropna()
Qin_data = jeme_basin.Qm3s.to_numpy() + 0.1
Qtime_data = jeme_basin.SOY.to_numpy()

#Import head water level (for comparison compare)
jeme_moulin = pd.read_csv('Field_Data/head_jeme.csv')
jeme_moulin = jeme_moulin.dropna()
h_real = jeme_moulin.head_bed.to_numpy()
t_real = jeme_moulin.soy.to_numpy()

#time parameters -- for 
time_start = Qtime_data[1200]
time_end = time_start + 50*secinday
timestep = 300 #seconds
time = TimeStamps(time_start,time_end,timestep)

meltwater_input = Qin_real(time, Qin_data, Qtime_data)
nsteps = len(meltwater_input)
baseflow = 0

# Qin_mean = 1
# dQ = 0.1 
# meltwater_input = Qin_sinusoidal(time,Qin_mean, dQ)


#%% Simulation with MouSh


moulin_evol = MoulinShape(                      
                    dz = dz,
                    z_elevations = z_elevations,
                    moulin_radii = moulin_radii,                  
                    temperature_profile = temperature_profile,                   
                    ice_thickness = ice_thickness,
                    initial_head = initial_head,
                    initial_subglacial_area = initial_subglacial_area,                           
                    regional_surface_slope = regional_surface_slope,
                    channel_length = channel_length,
                    creep_enhancement_factor = creep_enhancement_factor,                   
                    sigma = sigma,
                    tau_xy = tau_xy, 
                    friction_factor_OC = friction_factor_OC,
                    friction_factor_TM = friction_factor_TM,
                    friction_factor_SUB = friction_factor_SUB,
                    fraction_pd_melting = fraction_pd_melting,)

for idx,t in enumerate(time):
    moulin_evol.run1step(time,
                    timestep,
                    meltwater_input[idx],
                    subglacial_baseflow = baseflow, 
                    head_L = None,
                    overflow=True)   

                    
    
moulin_fix = MoulinShape(                      
                    dz = dz,
                    z_elevations = z_elevations,
                    moulin_radii = moulin_radii,                  
                    temperature_profile = temperature_profile,                   
                    ice_thickness = ice_thickness,
                    initial_head = initial_head,
                    initial_subglacial_area = initial_subglacial_area,                           
                    regional_surface_slope = regional_surface_slope)

for idx,t in enumerate(time):
    moulin_fix.run1step(t,
                        timestep, 
                        meltwater_input[idx],
                        creep=False,
                        elastic_deformation=False,
                        melt_below_head=False,
                        open_channel_melt=False,
                        potential_drop=False,
                        ice_motion=False,
                        refreezing=False,
                        overflow=True,
                        subglacial_baseflow = baseflow)



#%% Simulation with discretized moulin
discretized = run1Dsim(nsteps=nsteps, #number of timesteps
                       nx = 50, #number of nodes
                       dt = timestep,
                       L = channel_length,
                       Ztype = 'sqrt', 
                       Qin_type = 'custom_array',
                       Qsteady = False,
                       Qtime_data = time-time_start, #start time at zero
                       Qin_data = meltwater_input,
                       bedslope = regional_surface_slope,
                       A_R = np.pi*moulin_radii**2, # moulin cross-section area
                       D0 = 2 *np.pi*moulin_radii**2 / (np.pi*moulin_radii+ moulin_radii*2), # initial hydraulic diameter of the conduit
                       hin = initial_head, # initial upstream head
                       hout = 0)

#plt.plot(time/secinday,discretized['h'][:,0])
#%% Plot head for all simulations

head_evol = moulin_evol.dict['head']
head_fix = moulin_fix.dict['head']
head_discr = discretized['h'][:,0]

plt.figure()

plt.subplot(2,1,1)
plt.plot(time/secinday,meltwater_input)
plt.xlim([time_start/secinday,time_end/secinday])
plt.ylabel('Qin $m^3/s$')

plt.subplot(2,1,2)
plt.plot(t_real/secinday,h_real, color='black',label='field data')
plt.plot(time/secinday,head_evol,label='0D - evol')
plt.plot(time/secinday,head_fix,label='0D - fix')
plt.plot(time/secinday,head_discr,label='1D - fix')
plt.xlim([time_start/secinday,time_end/secinday])
plt.ylabel('Head (m)')
plt.xlabel('Day of year 2017')
plt.legend()

# add subplot for subglacial radius?
# add subplot for moulin radius at h?
















