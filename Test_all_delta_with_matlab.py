import function_moulin_model as fmm
import matplotlib.pyplot as plt
#import numpy as np

H = 1000
hw = 800
Mr_top=1
Mr_bottom=1

dz = 1
z = fmm.generate_grid_z(H,dz) 
Mr = fmm.initiate_moulin_radius(z,type='linear',Mr_top=Mr_top,Mr_bottom=Mr_bottom)

#%% dE

sigma_x = 0e3#-50e3 #(Units??) compressive
sigma_y = 50e3#-50e3 #(Units??) compressive
tau_xy = -50e3#100e3 #(Units??) shear opening

wet = fmm.locate_water(hw,z) 
Pw_z = fmm.calculate_water_pressure_at_depth(hw,z,wet) #water pressure at each depth
Pi_z = fmm.calculate_ice_pressure_at_depth(H,z) #ice pressure at each depth
sigma_z = fmm.calculate_sigma_z(Pw_z, Pi_z)

dE = fmm.calculate_elastic_deformation(Mr, sigma_z, sigma_x, sigma_y, tau_xy)
plt.figure()
plt.plot(dE,z,label='dE')
plt.legend()

#%% dC

dt = 300
E = 5
T_ice = fmm.interpolate_T_profile(z,temperature_profile=[273.15,273.15])
wet = fmm.locate_water(hw,z) 
Pw_z = fmm.calculate_water_pressure_at_depth(hw,z,wet) #water pressure at each depth
Pi_z = fmm.calculate_ice_pressure_at_depth(H,z) #ice pressure at each depth
sigma_z = fmm.calculate_sigma_z(Pw_z, Pi_z)
iceflow_param_glen = fmm.calculate_iceflow_law_parameter(T_ice,Pi_z)

dC = fmm.calculate_creep_moulin(Mr,dt,iceflow_param_glen,sigma_z,E)
plt.figure()
plt.plot(dC,z,label='dC')
plt.legend()
#%% dTM

Qout = 2;
friction_factor_TM = 0.1
dt = 300

wet = fmm.locate_water(hw,z) 
[Mcs, Mpr, Mdh, Mrh] = fmm.calculate_moulin_geometry(Mr, Mr)
uw_TM = fmm.calculate_water_velocity(Qout, Mcs)#,wet)
dL = fmm.calculate_dL(Mr, dz)
head_loss_dz_TM = fmm.calculate_moulin_head_loss(uw_TM, friction_factor_TM, dL, Mdh)
Pw_z = fmm.calculate_water_pressure_at_depth(hw,z,wet) #water pressure at each depth
Tmw = fmm.calculate_pressure_melting_temperature(Pw_z) 
T_ice = fmm.interpolate_T_profile(z,temperature_profile=[273.15,273.15])

dTM = fmm.calculate_melt_below_head(dL,head_loss_dz_TM, dt, Qout, Mpr, wet,include_ice_temperature=True,T_ice=T_ice,Tmw=Tmw)
plt.figure()
plt.plot(dTM,z,label='dTM')
plt.legend()
#%% dOC

Qin = 2
dt = 300
friction_factor_OC = 0.5

wet = fmm.locate_water(hw,z) 
dL = fmm.calculate_dL(Mr, dz)
[Mcs, Mpr, Mdh, Mrh] = fmm.calculate_moulin_geometry(Mr, Mr)
uw_OC = fmm.calculate_water_velocity(Qin, Mcs)#,wet)
head_loss_dz_OC = fmm.calculate_moulin_head_loss(uw_OC, friction_factor_OC, dL, Mdh)
T_ice = fmm.interpolate_T_profile(z,temperature_profile=[273.15,273.15])

dOC = fmm.calculate_melt_above_head_OC(Mr,Mpr,dL,dt,head_loss_dz_OC,Qin,wet,include_ice_temperature=True,T_ice=T_ice)
plt.figure()
plt.plot(dOC,z,label='dOC')
plt.legend()
#%% dF

t=300
dt = 300
wet = fmm.locate_water(hw,z) 
T_ice = fmm.interpolate_T_profile(z,temperature_profile=[260,260])

dF = fmm.calculate_refreezing(T_ice,t,dt,wet)
dF_test = fmm.calculate_refreezing_test(T_ice,t,dt,wet)
plt.figure()
plt.plot(dF,z,label='dF')
plt.legend()
plt.figure()
plt.plot(dF_test,z,label='dF test')
plt.legend()

