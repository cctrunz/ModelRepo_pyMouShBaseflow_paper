''' 
Moulin Physical Model from LC Andrews and K Poinar
Translated by Celia Trunz

Component provenance:
    - Subglacial channel: Schoof2010
    - 
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import function_moulin_model as fmm

import pandas as pd




secinday=24*3600


'''Default parameter'''

H = 1000 #(m) Ice thickness
L = 10000
Mr_top=0.1
Mr_bottom=0.1
regional_surface_slope =0.0# alpha in matlab
E = 5 
sigma_x = 0e3#-50e3 #(Units??) compressive
sigma_y = 50e3#-50e3 #(Units??) compressive
tau_xy = -50e3#100e3 #(Units??) shear opening

tmax_in_day = 10#600 /3600/24 #(days) Maximum time to run
dt = 300 #(s) timestep
time = fmm.generate_time(dt,tmax_in_day)

dz = 1 #(m) 
z = fmm.generate_grid_z(H,dz)

Qin_mean = 1 #m3/s
period = 24*3600
dQ = 0.2
Q_type = 'sinusoidal_celia'
Qin = fmm.set_Qin(time,type=Q_type,Qin_mean=Qin_mean,dQ=dQ,period=period)

T_ice = fmm.interpolate_T_profile(z,temperature_profile=[273.15,273.15])



fraction_pd_melting = 0 #This is fraction of energy from potential drop that is used for melting
friction_factor_TM = 0.1
friction_factor_OC = 0.5


Pi_H = fmm.calculate_ice_pressure_at_depth(H,0) 
Pi_z = fmm.calculate_ice_pressure_at_depth(H,z) 
iceflow_param_glen = fmm.calculate_iceflow_law_parameter(T_ice,Pi_z) 

#Initialize

hw = np.zeros(len(time))
hw[0] = H #(m)Initial water level
SCs = 1 

Mr_major = fmm.initiate_moulin_radius(z,type='linear',Mr_top=Mr_top,Mr_bottom=Mr_bottom)
Mr_minor = Mr_major
[Mx_upstream, Mx_downstream]= fmm.initiate_moulin_wall_position(Mr_major,Mr_minor) #for position of the wall relative to coordinate system

dGlen = 0
dGlen_cumulative = 0
Vadd_C = 0
Vadd_E = 0
Vadd_TM = 0



for idx, t in enumerate(time):
    
      
    '''Calculate moulin geometry, relative to water level'''    
    [Mcs, Mpr, Mdh, Mrh] = fmm.calculate_moulin_geometry(Mr_major, Mr_minor)
    
    '''Calculate water level''' 
    Qin_compensated = Qin[idx]#+ Vadd_E + Vadd_C
    sol = solve_ivp(fmm.calculate_h_S_schoof,
                    [0, dt], #initial time and end time. We solve for one timestep.
                    [hw[idx-1],SCs], #initial head and channel cross-section area. Uses values in previous timestep.
                    args=(Mcs,z,Pi_H,L,Qin_compensated,H,False), #args() only works with scipy>=1.4. if below, then error message: missing 5 arguments
                    method = 'LSODA' #solver method
                    # atol = 1e-6, #tolerance. can make it faster
                    # rtol = 1e-3,
                    #max_step = 10 #change if resolution is not enough
                    )
    hw[idx] = sol.y[0][-1]  #(m) moulin water level
    SCs = sol.y[1][-1] #(m) Channel cross-section
    Qout = fmm.calculate_Qout(SCs,hw[idx],L)

    '''update parameters (in loop)'''
    wet = fmm.locate_water(hw[idx],z) 
    Pw_z = fmm.calculate_water_pressure_at_depth(hw[idx],z,wet) #water pressure at each depth
    Tmw = fmm.calculate_pressure_melting_temperature(Pw_z)   
    sigma_z = fmm.calculate_sigma_z(Pw_z, Pi_z)
    uw_TM = fmm.calculate_water_velocity(Qout, Mcs)
    uw_OC = fmm.calculate_water_velocity(Qin[idx], Mcs)
    dL_upstream = fmm.calculate_dL(Mx_upstream, dz)
    head_loss_dz_TM = fmm.calculate_moulin_head_loss(uw_TM, friction_factor_TM, dL_upstream, Mdh)
    head_loss_dz_OC = fmm.calculate_moulin_head_loss(uw_OC, friction_factor_OC, dL_upstream, Mdh)
    
    
    '''Calculate moulin changes for each component'''
    #Creep Deformation
    dC_major = fmm.calculate_creep_moulin(Mr_major,dt,iceflow_param_glen,sigma_z,E)
    dC_minor = fmm.calculate_creep_moulin(Mr_major,dt,iceflow_param_glen,sigma_z,E)
    Vadd_C = fmm.calculate_Q_stress_wall(dC_major,dC_minor,Mr_major,Mr_minor,z,wet,dt)    
    #Turbulent melting
    dTM = fmm.calculate_melt_below_head(dL_upstream,head_loss_dz_TM, dt, Qout, Mpr, wet,include_ice_temperature=True,T_ice=T_ice,Tmw=Tmw)
    vadd_TM = fmm.calculate_Q_melted_wall(dTM, z, Mpr, dt)    
    #Refreezing
        
    #Open channel melting
    dOC = fmm.calculate_melt_above_head_OC(Mr_major,Mpr,dL_upstream,dt,head_loss_dz_OC,Qin[idx],wet,include_ice_temperature=True,T_ice=T_ice)
    vadd_OC = fmm.calculate_Q_melted_wall(dOC, z, Mpr/2, dt) 
    dPD = fmm.calculate_melt_above_head_PD(Mr_major, Qin[idx], dt, Mpr, wet, fraction_pd_melting)   
    vadd_PD = fmm.calculate_Q_melted_wall(dPD, z, Mpr/2, dt) 
    #Elastic deformation
    dE_major = fmm.calculate_elastic_deformation(Mr_major, sigma_z, sigma_x, sigma_y, tau_xy)
    dE_minor = fmm.calculate_elastic_deformation(Mr_minor, sigma_z, sigma_x, sigma_y, tau_xy)
    Vadd_E = fmm.calculate_Q_stress_wall(dE_major,dE_minor,Mr_major,Mr_minor,z,wet,dt)
    #Asymmetric deformation due to Glen's Flow Law
    dGlen = fmm.calculate_iceflow_moulin(iceflow_param_glen, regional_surface_slope, H, z, dt)
    dGlen_cumulative = fmm.calculate_cumulative_dGlen(dGlen, dGlen_cumulative)
    
    
    '''Update moulin radii'''   
    dr_major = fmm.calculate_dradius(dTM=dTM)#dC=dC_major)#dE=dE_major)#, dC=dC_major, dTM=dTM)#, dOC=dOC)#dPD=dPD
    dr_minor = fmm.calculate_dradius(dTM=dTM)#dC=dC_major)#dE=dE_minor)#, dC=dC_major, dTM=dTM)
    Mr_major = fmm.update_moulin_radius( Mr_major,dr_major )
    Mr_minor = fmm.update_moulin_radius( Mr_minor,dr_minor )
    [Mx_upstream, Mx_downstream] = fmm.update_moulin_wall_position(Mx_upstream, Mx_downstream, dr_major,dr_minor, dGlen, dGlen_cumulative)



#%%
plt.figure()
plt.plot(Mx_upstream,z)
plt.plot(Mx_downstream,z)
plt.xlim([-3,3])

#%%
plt.figure()
plt.plot(time,hw)
