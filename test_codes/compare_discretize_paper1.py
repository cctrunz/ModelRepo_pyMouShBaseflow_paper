# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 15:54:52 2021

@author: celia
"""

# def Qin_sinusoidal(time,Qin_mean, dQ, shift=0):
#     """
#     dQ - amplitude of meltwater oscillation
#     Qin_mean
#     """
#     # Qin_mean=3, dQ=0.5, period=24*3600
#     return dQ*np.sin(2*np.pi*time/secinday+shift) + Qin_mean

## Compare Moush (0D), fixed cylinder (0D), and fixed cylinder (1D)
import numpy as np
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import seaborn as sns


#for the 1D simulation. This requires to have Matt's modules installed
from conduits1D_landlab_matt_02112021 import run1Dsim,plot_3panels_singlestep, plot_3panels_movie, plot_2panel_overtime_multiposition #for the 1D simulation. This requires to have Matt's modules installed

from pyMouSh import MoulinShape, TimeStamps, Qin_sinusoidal, Qin_real

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

#Glacier parameters
moulin_radii = 1. #m (based out of the highest radius when running moush for JEME)
ice_thickness = 1000 #m
initial_head = 900 #m
initial_subglacial_area = 1 #m
channel_length = 30e3 #m 
friction_factor_SUB = 0.1 # its the same value in pressurized_flow in landlab
regional_surface_slope = 0

#time parameters -- for 
time_start = 0
time_end = 20*secinday
timestep = 5*60 #300 #seconds
time = TimeStamps(time_start,time_end,timestep)

Qin_mean = 3
dQ = 0.3
meltwater_input = Qin_sinusoidal(time,Qin_mean, dQ, shift=0) #Qin_real(time, Qin_data, Qtime_data)

nsteps = len(meltwater_input)


#%% Simulation with MouSh
    
moulin_fix = MoulinShape(                      
                    dz = dz,
                    z_elevations = z_elevations,
                    moulin_radii = moulin_radii,                  
                    temperature_profile = temperature_profile,                   
                    ice_thickness = ice_thickness,
                    initial_head = initial_head,
                    channel_length = channel_length,
                    friction_factor_SUB = friction_factor_SUB,
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
                        subglacial_baseflow = 0)

# picklefile = open('moulin_fix', 'wb')
# pickle.dump(moulin_fix, picklefile)
# picklefile.close()

#%% Simulation with discretized moulin
radius = moulin_radii

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

# picklefile = open('discretized', 'wb')
# pickle.dump(discretized, picklefile)
# picklefile.close()

#plt.plot(time/secinday,discretized['h'][:,0])
#%% Plot head for all simulations
#xlim = [195,200]
#xlim = [180,250]
#☻xlim = [0,20]

head_fix = moulin_fix.dict['head']
#head_discr = discretized['h'][:,0]

r_fix = np.array(moulin_fix.dict['subglacial_radius'])

dh_fix = 2 *np.pi*r_fix**2 / (np.pi*r_fix+ r_fix*2)

fig,ax = plt.subplots(3,1, sharex=True)
#plt.figure()

ax[0].plot(time/secinday,meltwater_input, color='grey')
#plt.xlim(xlim)
ax[0].set_ylabel('Q$_{in}$ (m$^3$/s)')
ax[0].tick_params(labelbottom=False)
ax[0].set_yticks([2,3,4])


ax[1].plot(time/secinday,head_fix, linestyle='--',label='0D (moulin)')
ax[1].plot(time/secinday,discretized['h'][:,0],label='1D (moulin)')
ax[1].plot(time/secinday,discretized['h'][:,-2],label='1D (margin)')
ax[1].set_xlim([time_start/secinday,time_end/secinday])
ax[1].set_ylabel('h (m)')
ax[1].legend(loc=4, prop={'size': 6}, bbox_to_anchor=(1,-0.2))
ax[1].tick_params(labelbottom=False)


ax[2].plot(time/secinday, dh_fix,linestyle='--',label='0D (moulin)')
ax[2].plot(time/secinday,discretized['d'][:,0],label='1D (moulin)')
ax[2].plot(time/secinday,discretized['d'][:,-2],label='1D (margin)')
ax[2].set_xlim([time_start/secinday,time_end/secinday])
ax[2].set_yticks([1,1.2,1.4])
ax[2].set_xticks([0,2,4,6,8,10,12,14,16,18,20])
ax[2].set_ylabel('HD (m)')
ax[2].set_xlabel('Days')

ax[0].text(0.3,3.8,'(a)',fontsize=8)
ax[1].text(0.3,1100,'(b)',fontsize=8)
ax[2].text(0.3,1.3,'(c)',fontsize=8)

sns.despine(trim=True)

plt.savefig('compare_0D1D.pdf')
plt.savefig('compare_0D1D.png')

# add subplot for subglacial radius?
# add subplot for moulin radius at h?
















