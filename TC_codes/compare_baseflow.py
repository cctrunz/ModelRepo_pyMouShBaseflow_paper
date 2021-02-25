
from pyMouSh.pyMouSh import MoulinShape, TimeStamps, Qin_constant, Qin_sinusoidal, Qin_real, calculate_h_S_schoof
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle

#%%
secinday = 24*3600
ZERO_KELVIN = 273.15
#paramters
                                
dz = 1
z_elevations = None
moulin_radii = 0.3
temperature_profile = np.array([ZERO_KELVIN, ZERO_KELVIN])
ice_thickness = 500
initial_head = 500
initial_subglacial_area = 1 #back to 1!!!      
regional_surface_slope = 0
channel_length = 25e3
creep_enhancement_factor = 3
sigma = (0, 0)  
tau_xy = 0  
friction_factor_OC = 1
friction_factor_TM = 0.5
friction_factor_SUB = 0.1
fraction_pd_melting = 0.2


jeme_basin = pd.read_csv('MouSh_baseflow_code_figures/Field_Data/surface_melt_jeme.csv')
jeme_basin = jeme_basin.dropna()
Qin_data = jeme_basin.Qm3s.to_numpy() +0.1
Qtime_data = jeme_basin.SOY.to_numpy() + 3*3600 #add routing delay -- Jess found 2h not 4. investigate why


time_start = Qtime_data[int(2*secinday/900)]  
time_end = time_start + 50*secinday
timestep = 30*60#300 #seconds to reduce number of saved data
time = TimeStamps(time_start,time_end,timestep)

meltwater_input = Qin_real(time, Qin_data, Qtime_data)

 

#%%

bf0_fix = MoulinShape(                      
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
                    fraction_pd_melting = fraction_pd_melting)

for idx,t in enumerate(time) :
    bf0_fix.run1step(time,
                    timestep,
                    meltwater_input[idx],
                    subglacial_baseflow = 0, 
                    head_L = None )

print('bf0 finished, starting pickle')

picklefile = open('bf0_fix', 'wb')
pickle.dump(bf0_fix, picklefile)
picklefile.close()

print('bf0 finished')

del bf0_fix

#%%

bf3_fix = MoulinShape(                      
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
                    fraction_pd_melting = fraction_pd_melting)

for idx,t in enumerate(time) :
    bf3_fix.run1step(time,
                    timestep,
                    meltwater_input[idx],
                    subglacial_baseflow = 3, 
                    head_L = None )

print('bf3_fix finished, starting pickle')

picklefile = open('bf3_fix', 'wb')
pickle.dump(bf3_fix, picklefile)
picklefile.close()

print('bf3_fix finished')

del bf3_fix

#%%

bf1_fix = MoulinShape(                      
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
                    fraction_pd_melting = fraction_pd_melting)

for idx,t in enumerate(time) :
    bf1_fix.run1step(time,
                    timestep,
                    meltwater_input[idx],
                    subglacial_baseflow = 1, 
                    head_L = None )

print('bf1_fix finished, starting pickle')

picklefile = open('bf1_fix', 'wb')
pickle.dump(bf1_fix, picklefile)
picklefile.close()

print('bf1_fix finished')

del bf1_fix
#%%

Qin_mean = 1
dQ = 0.1 
shift1 = 0 * secinday # sin is already shifted from input
baseflow_shift = Qin_sinusoidal(time,Qin_mean, dQ, shift=shift1)

bf1_osc_shift = MoulinShape(                      
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
                    fraction_pd_melting = fraction_pd_melting)

for idx,t in enumerate(time) :
    bf1_osc_shift.run1step(time,
                    timestep,
                    meltwater_input[idx],
                    subglacial_baseflow = baseflow_shift[idx], 
                    head_L = None )

print('bf1_osc_shift finished, starting pickle')

picklefile = open('bf1_osc_shift', 'wb')
pickle.dump(bf1_osc_shift, picklefile)
picklefile.close()

print('bf1_osc_shift finished')

del bf1_osc_shift

#%%

Qin_mean = 1
dQ = 0.1 
shift2 = 0.42 * secinday
baseflow_osc= Qin_sinusoidal(time,Qin_mean, dQ, shift=shift2)

bf1_osc = MoulinShape(                      
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
                    fraction_pd_melting = fraction_pd_melting)

for idx,t in enumerate(time) :
    bf1_osc.run1step(time,
                    timestep,
                    meltwater_input[idx],
                    subglacial_baseflow = baseflow_osc[idx], 
                    head_L = None )

print('bf1_osc finished, starting pickle')

picklefile = open('bf1_osc', 'wb')
pickle.dump(bf1_osc, picklefile)
picklefile.close()

print('bf1_osc finished')

del bf1_osc


#%%


baseflow_melt_shift = Qin_real(time + 12*3600, Qin_data, Qtime_data)

bf_melt_shift = MoulinShape(                      
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
                    fraction_pd_melting = fraction_pd_melting)

for idx,t in enumerate(time) :
    bf_melt_shift.run1step(time,
                    timestep,
                    meltwater_input[idx],
                    subglacial_baseflow = baseflow_melt_shift[idx], 
                    head_L = None )

print('bf_melt_shift finished, starting pickle')

picklefile = open('bf_melt_shift', 'wb')
pickle.dump(bf_melt_shift, picklefile)
picklefile.close()

print('bf_melt_shift finished')

del bf_melt_shift

