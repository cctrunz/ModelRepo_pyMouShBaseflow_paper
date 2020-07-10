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
import plot_moulin_model as pmm

'''Activate or deactivate False and True statements'''

include_ice_temperature = False #%true means that the change in the ice temperature is included in...
#%the calculated change in moulin radius. If false, it makes the implicit
#%assumption that the ice temperature and water temperature are both at the pressure melting temperature.

#do the same for Bathurst....




'''Default parameter'''
#R0 = 2 #(m) Initial moulin radius
H = 1000 #(m) Ice thickness
L = 10000 #(m) Subglacial channel length
tmax_in_day = 5 #(days) Maximum time to run
dt = 300 #(s) timestep
Mr_minor_initial = 1 #(m)
Mr_major_initial = 1 #(m)
Mr_minimum = 1e-9 #(m)
xmax    = 30 # 80 #(m) how far away from moulin to use as infinity
hw = H #(m)Initial water level
SCs = 1.5 #(m) Initial subglacial channel cross-section area
Qin_mean = 2.5 #m3/s
dQ = 1
E = 5 #Enhancement factor for the ice creep.
regional_surface_slope = 0.01#alpha in matlab %regional surface slope (unitless), for use in Glen's Flow Law
#Q_type = 'constant' #
Q_type = 'sinusoidal_celia'
T_type = 'Cool' 

#Assign elastic deformation parameters
#stress={} #initiate dictionnary for elastic deformation parameters
sigma_x = -50e3 #(Units??) compressive
sigma_y = -50e3 #(Units??) compressive
tau_xy = 100e3 #(Units??) shear opening

#Turbulent melting parameters
relative_roughness = 1 #increasing this value increases the amount of melting due to turbulence. (comment from Matlab code)
relative_roughness_OC = 1e-9 #1e-12;  % This one modifies the melt from open channel flow.

#initialize
dGlen = 0
dGlen_cumulative = 0

'''set and initialize parameters and dimensions'''
#---set model grid
[z, nz, dz] = fmm.generate_grid_z(H) #default is 1m spacing # should this mimic the x grid??
#---set moulin geom 
#[Mr_major, Mr_minor]= fmm.initiate_moulin_radius(z,Mr_major_initial,Mr_minor_initial) #for moulin dimension
[Mx_upstream, Mx_downstream, Mr_major, Mr_minor]= fmm.initiate_moulin_wall_position(Mr_major_initial, Mr_minor_initial,z) #for position of the wall relative to coordinate system
[x, nx, dx] = fmm.generate_grid_x(dt, xmax) #generate grid around the moulin for temperature discretization
#---set duration of the model run'
[time,tmax_in_second] = fmm.generate_time(dt,tmax_in_day)
#---set ice temperature
[T_far,T] = fmm.set_ice_temperature(x,z,type=T_type)
#specific for turbulence
if include_ice_temperature:
    Ti = T_far
else:
    Ti = np.nan
    
'''Calculate or import meltwater input'''
#sinusoidal_celia
#constant
Qin = fmm.set_Qin(time,type=Q_type,Qin_mean=Qin_mean,dQ=dQ)

'''Calculate parameters (out of the loop)'''
#---calculate total ice pressure
Pi_H = fmm.calculate_ice_pressure_at_depth(H,0) #total ice pressure (ice pressure at the bottom)
Pi_z = fmm.calculate_ice_pressure_at_depth(H,z) #ice pressure at each depth
T_mean = fmm.calculate_mean_ice_temperature(T)
iceflow_param_glen = fmm.calculate_iceflow_law_parameter(T_mean,Pi_z) #(units?) called A in matlab's code
#stress_cryo = -Pi_z # Ice hydrostatic stress (INWARD: Negative)

#Initiate results dictionnary
results = fmm.initiate_results_dictionnary(time,z)
#Save parameter in dictionnary


#makes an artificial cone
# Mr_major=Mr_major*np.linspace(1,2,len(z))
# Mr_minor=Mr_minor*np.linspace(1,2,len(z))



for idx, t in enumerate(time):
    

    
    '''Calculate moulin geometry, relative to water level'''
    
    [Mcs, Mpr, Mdh, Mrh,Diameter] = fmm.calculate_moulin_geometry(
        Mx_upstream, Mx_downstream, Mr_major, Mr_minor)
    
    '''Calculate water level''' 
    sol = solve_ivp(fmm.calculate_h_S_schoof,
                    [0, dt], #initial time and end time. We solve for one timestep.
                    [hw,SCs], #initial head and channel cross-section area. Uses values in previous timestep.
                    args=(Mcs,z,Pi_H,L,Qin[idx],H,False), #args() only works with scipy>=1.4. if below, then error message: missing 5 arguments
                    method = 'LSODA' #solver method
                    # atol = 1e-6, #tolerance. can make it faster
                    # rtol = 1e-3,
                    #max_step = 10 #change if resolution is not enough
                    )
    hw = sol.y[0][-1]  #(m) moulin water level
    SCs = sol.y[1][-1] #(m) Channel cross-section
    Qout = fmm.calculate_Qout(SCs,hw,L)


    '''Calculate parameters (in loop)'''


    wet = fmm.locate_water(hw,z) 
    Pw_z = fmm.calculate_water_pressure_at_depth(hw,z,wet) #water pressure at each depth
    Tmw = fmm.calculate_pressure_melting_temperature(Pw_z)   
    #stress_hydro = Pw_z # Water hydrostatic stress (OUTWARD: Positive)'
    Pwi_z = Pw_z - Pi_z # water pressure open,  
    uw = fmm.calculate_water_velocity(Qout, Mcs, wet)
    fR_bathurst = fmm.calculate_bathurst_friction_factor(Mrh, relative_roughness)
    fR_colebrook_white = fmm.calculate_colebrook_white_friction_factor(Mdh, relative_roughness)
    
    
    
    '''Calculate moulin changes for each component'''
    #Creep Deformation

    [dC_major,dC_minor] = fmm.calculate_creep_moulin(Mr_major,Mr_minor,dt,iceflow_param_glen,Pwi_z,E)

    
    #Turbulent melting
    dTM = fmm.calculate_turbulent_melting_moulin(
        Mx_upstream, Mx_downstream, fR_bathurst, uw, Tmw, Ti, z, dz, dt, Qout, Mpr, Mdh, include_ice_temperature)

    vadd_turb = fmm.calculate_volume_melted_wall(dTM, z, Mpr, dt)
    
    #Refreezing
    
    
    #Open channel melting
    dOC = fmm.calculate_melt_above_head(Mr_major, Qin[idx], dt, Mpr, wet, method='potential_drop')
    
    #Elastic deformation
    [dE_major,dE_minor] = fmm.calculate_elastic_deformation(Mr_major, Mr_minor, Pwi_z, sigma_x, sigma_y, tau_xy)
    #Asymmetric deformation due to Glen's Flow Law
    dGlen = fmm.calculate_iceflow_moulin(Pi_z, iceflow_param_glen, regional_surface_slope, H, z, dt)
    dGlen_cumulative = fmm.calculate_cumulative_dGlen(dGlen, dGlen_cumulative)
    
    
    '''Calculate new moulin radii'''
 
    
    [Mx_upstream, Mx_downstream, Mr_major, Mr_minor] = fmm.calculate_new_moulin_wall_position( Mx_upstream, Mx_downstream, dGlen_cumulative,
                                                                    dC = [dC_major, dC_minor],
                                                                    dTM = dTM,
                                                                    dGlen = dGlen,
                                                                    dE = [dE_major,dE_minor],
                                                                    dOC = dOC
                                                                    )
    
    '''Save values'''
    results['Mx_upstream'][idx] = Mx_upstream
    results['Mx_downstream'][idx] = Mx_downstream
    results['Mr_major'][idx] = Mr_major
    results['Mr_minor'][idx] = Mr_minor
    results['Diameter'][idx] = Diameter
    results['dC_major'][idx] = dC_major
    results['dC_minor'][idx] = dC_minor
    results['dTM'][idx] = dTM
    results['dOC'][idx] = dOC
    # results['dTM_major'][idx] = dTM_major
    # results['dTM_minor'][idx] = dTM_minor
    results['dGlen'][idx] =  dGlen
    results['dGlen_cumulative'][idx] =  dGlen_cumulative
    results['dE_major'][idx] = dE_major
    results['dE_minor'][idx] = dE_minor
    results['Mcs'][idx] =  Mcs
    results['Mpr'][idx] = Mpr
    results['Mdh'][idx] = Mdh
    results['Mrh'][idx] = Mrh
    results['Pi_z'][idx] = Pi_z
    results['Pw_z'][idx] = Pw_z
    results['wet'][idx] = wet
    results['Tmw'][idx] = Tmw
    results['Pwi_z'][idx] = Pwi_z
    results['uw'][idx] = uw
    results['fR_bathurst'][idx] = fR_bathurst
    results['fR_colebrook_white'][idx] = fR_colebrook_white
    
    results['hw'][idx] = hw
    results['SCs'][idx] = SCs
    results['Qout'][idx] = Qout
    
    
    
'''Plot'''
pmm.plot_geom(results, time, z, set_xlim_moulin=False,
              ax2_varname=['Mr_major','Mr_minor'],
              ax3_varname=['dC_major','dC_minor'],
              ax4_varname='dTM',
              ax5_varname=['dE_major','dE_minor'], 
              ax6_varname='dOC' )
#pmm.plot_geom(results,ax2_varname=['dC_major','dC_minor'],ax3_varname='Pw_z',ax4_varname='Pi_z',ax5_varname='Pwi_z'  )
#pmm.plot_2Darray(results,results['dTM'])

#pmm.plot_2Darray_with_1Darray(results,results['Mr_major'],results['hw'])

pmm.plot_1Darray_timeserie(results, time, 'hw')
pmm.plot_1Darray_timeserie(results, time, 'SCs')
plt.figure()
plt.plot(time,Qin)



# colors = [plt.cm.rainbow(i) for i in np.linspace(0, 1, len(results['time']))] 
# plt.figure()
# for i in np.arange(len(time)):
#     #plt.plot(results['Mr_major'][i],results['z']) 
#     plt.plot(results['Pi_z'][i],results['z'],color=colors[i]) 
#     plt.plot(results['Pw_z'][i],results['z'],color=colors[i]) 
    
# plt.figure()    
# d_Mr = results['Mr_major']-results['Mr_major']   
# for i in np.arange(0,len(results['time']),100):   
#     plt.plot(d_Mr[i],z)
    
# plt.figure()    
# for i in np.arange(0,len(results['time']),100):   
#     plt.plot(results['Mr_minor'][i],z)   


