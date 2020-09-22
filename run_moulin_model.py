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
import pandas as pd
#import comprehensive_plot
#import plot_pretty_moulin
#import plot_deltas
import plot_all_in_one
#import plot_pretty_moulin_with_deltas



secinday=24*3600


'''Default parameter'''
#R0 = 2 #(m) Initial moulin radius
H = 1000 #(m) Ice thickness
#-- can be fixed or can be calculated for an idealized profile
L = fmm.calculate_L(H) #L = 10000 #(m) Subglacial channel length 
tmax_in_day = 10 #(days) Maximum time to run
dt = 300 #(s) timestep
mts_to_cmh = 100*60*60/dt #m per timestep to mm/h : change units
Mr_minor_initial = 1 #(m)
Mr_major_initial = 1 #(m)
Mr_minimum = 1e-9 #(m)
xmax    = 30 # 80 #(m) how far away from moulin to use as infinity
hw = H #(m)Initial water level
SCs = 1.5 #(m) Initial subglacial channel croohhhss-section area
Qin_mean = 1 #m3/s
dQ = 0.1
E = 3 #Enhancement factor for the ice creep.
regional_surface_slope = 0.02#alpha in matlab %regional surface slope (unitless), for use in Glen's Flow Law
#Q_type = 'constant' #
Q_type = 'sinusoidal_celia'
#import some temperature profiles from the litterature											
temperature_profil_litterature = pd.read_excel('FieldData/Temp_Profil_Litterature.xlsx',sheet_name=None) #import all the sheet at the same time
print(temperature_profil_litterature.keys()) #list the temperature profiles options

Temperature_profile = temperature_profil_litterature['Lüthi15_FOXX1'].temperature #[273.15,273.15]
#reverse data to have the 0 m elevation z match the 0 row
Temperature_profile = Temperature_profile.iloc[::-1]
#Ryser_foxx is Luthi15_FOXX1 in matlab
#Ryser_GULL is Luthi15_GULL
# REF: 
    #Iken, A., Echelmeyer, K., Harrison, W. D. & Funk, M. Mechanisms of fast flow in Jakobshavns Isbræ, West Greenland: Part I. Measurements of temperature and water level in deep boreholes. Journal of Glaciology 39, 15–25 (1993).
    #Lüthi, M. P. et al. Heat sources within the Greenland Ice Sheet: dissipation, temperate paleo-firn and cryo-hydrologic warming. The Cryosphere 9, 245–253 (2015).
# or create array: temperate profile: [273.15,273.15], linear cold profile: [273.15,265]


#Assign elastic deformation parameters
#stress={} #initiate dictionnary for elastic deformation parameters
sigma_x = 0e3#-50e3 #(Units??) compressive
sigma_y = 50e3#-50e3 #(Units??) compressive
tau_xy = -50e3#100e3 #(Units??) shear opening

#Turbulent melting parameters
relative_roughness = 10 #increasing this value increases the amount of melting due to turbulence. (comment from Matlab code)
relative_roughness_OC = 1e-7 #1e-12;  % This one modifies the melt from open channel flow.

#initialize
dGlen = 0
dGlen_cumulative = 0
Vadd_C = 0
Vadd_E = 0
Vadd_TM = 0

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
[T_far,T] = fmm.set_ice_temperature(x,z,Temperature=Temperature_profile)
#specific for turbulence
include_ice_temperature = False #%true means that the change in the ice temperature is included in...
#%the calculated change in moulin radius. If false, it makes the implicit
#%assumption that the ice temperature and water temperature are both at the pressure melting temperature.
#do the same for Bathurst....
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

#makes an artificial cone
# Mr_major=Mr_major*np.linspace(1,2,len(z))
# Mr_minor=Mr_minor*np.linspace(1,2,len(z))

idx_plot = 0
i_save = 0




for idx, t in enumerate(time):
      
    '''Calculate moulin geometry, relative to water level'''    
    [Mcs, Mpr, Mdh, Mrh,Diameter] = fmm.calculate_moulin_geometry(Mx_upstream, Mx_downstream, Mr_major, Mr_minor)
    
    '''Calculate water level''' 
    Qin_compensated = Qin[idx]+ Vadd_E + Vadd_C
    sol = solve_ivp(fmm.calculate_h_S_schoof,
                    [0, dt], #initial time and end time. We solve for one timestep.
                    [hw,SCs], #initial head and channel cross-section area. Uses values in previous timestep.
                    args=(Mcs,z,Pi_H,L,Qin_compensated,H,False), #args() only works with scipy>=1.4. if below, then error message: missing 5 arguments
                    method = 'LSODA' #solver method
                    # atol = 1e-6, #tolerance. can make it faster
                    # rtol = 1e-3,
                    #max_step = 10 #change if resolution is not enough
                    )
    hw = sol.y[0][-1]  #(m) moulin water level
    SCs = sol.y[1][-1] #(m) Channel cross-section
    Qout = fmm.calculate_Qout(SCs,hw,L)

    '''update parameters (in loop)'''
    wet = fmm.locate_water(hw,z) 
    Pw_z = fmm.calculate_water_pressure_at_depth(hw,z,wet) #water pressure at each depth
    Tmw = fmm.calculate_pressure_melting_temperature(Pw_z)   
    #stress_hydro = Pw_z # Water hydrostatic stress (OUTWARD: Positive)'
    sigma_z = fmm.calculate_sigma_z(Pw_z, Pi_z)
    uw = fmm.calculate_water_velocity(Qout, Mcs, wet)
    fR_bathurst = fmm.calculate_bathurst_friction_factor(Mrh, relative_roughness)
    fR_colebrook_white = fmm.calculate_colebrook_white_friction_factor(Mdh, relative_roughness)   
    
    '''Calculate moulin changes for each component'''
    #Creep Deformation
    [dC_major,dC_minor] = fmm.calculate_creep_moulin(Mr_major,Mr_minor,dt,iceflow_param_glen,sigma_z,E)
    Vadd_C = fmm.calculate_Q_stress_wall(dC_major,dC_minor,Mr_major,Mr_minor,z,wet,dt)    
    #Turbulent melting
    dTM = fmm.calculate_melt_below_head(Mx_upstream, Mx_downstream, fR_bathurst, uw, Tmw, Ti, z, dz, dt, Qout, Mpr, Mdh, include_ice_temperature,wet)
    vadd_TM = fmm.calculate_Q_melted_wall_below_head(dTM, z, Mpr, dt)    
    #Refreezing
        
    #Open channel melting
    dPD = fmm.calculate_melt_above_head_PD(Mr_major, Qin[idx], dt, Mpr, wet)    
    #Elastic deformation
    [dE_major,dE_minor] = fmm.calculate_elastic_deformation(Mr_major, Mr_minor, sigma_z, sigma_x, sigma_y, tau_xy)
    Vadd_E = fmm.calculate_Q_stress_wall(dE_major,dE_minor,Mr_major,Mr_minor,z,wet,dt)
    #Asymmetric deformation due to Glen's Flow Law
    dGlen = fmm.calculate_iceflow_moulin(Pi_z, iceflow_param_glen, regional_surface_slope, H, z, dt)
    dGlen_cumulative = fmm.calculate_cumulative_dGlen(dGlen, dGlen_cumulative)
    
    
    
    
    '''Update moulin radii'''   
    [dr_major,dr_minor] = fmm.calculate_dradius(dC=[dC_major, dC_minor], dTM=dTM, dE=[dE_major,dE_minor], dPD=dPD)   
    [Mx_upstream, Mx_downstream, Mr_major, Mr_minor] = fmm.calculate_new_moulin_wall_position(Mx_upstream, Mx_downstream,Mr_major, Mr_minor, dr_major,dr_minor, dGlen, dGlen_cumulative)
    
    # if idx_plot == idx:
    #     i_save = i_save+1
    #     idx_plot = idx_plot+20
    #     plot_all_in_one.live_plot(Mx_upstream,Mx_downstream,dTM,dPD,dC_major,dC_minor,dE_major,dE_minor,dr_major,dr_minor,dGlen,t,hw,SCs,z,Qin,Qout,idx,wet,mts_to_cmh,time,H,results,T_far)
        
        #so for instance if “nn” is your run count
        
        # if mod(nn,10)
        # ((Plot  something))
        # end
        # OOPS
        # If ~mod(nn,10)
        # end
        
        #comprehensive_plot.live_plot(Mx_upstream,Mx_downstream,dTM,dPD,dC_major,dC_minor,dE_major,dE_minor,dr_major,dr_minor,dGlen,t,hw,SCs,z,Qin,Qout,idx,wet,mts_to_cmh,time,H)
        #plot_pretty_moulin.live_plot(hw,Mx_upstream,Mx_downstream,z)
        #plot_deltas.live_plot(dTM,dPD,dC_major,dC_minor,dE_major,dE_minor,dr_major,dr_minor,dGlen,z,wet,mts_to_cmh)
    #plot_pretty_moulin_with_deltas.live_plot(Mx_upstream,Mx_downstream,dTM,dPD,dC_major,dC_minor,dE_major,dE_minor,dr_major,dr_minor,dGlen,t,hw,z,idx,wet,mts_to_cmh,time,H)
        #plt.savefig('figure_movie/all_in_one_test/all_in_one_test%s.png'%i_save)
        #plt.pause(0.001)

    '''Save values'''
    results['Mx_upstream'][idx] = Mx_upstream
    results['Mx_downstream'][idx] = Mx_downstream
    results['Mr_major'][idx] = Mr_major
    results['Mr_minor'][idx] = Mr_minor
    results['Diameter'][idx] = Diameter
    results['dC_major'][idx] = dC_major
    results['dC_minor'][idx] = dC_minor
    results['dTM'][idx] = dTM
    #results['dOC'][idx] = dOC
    results['dPD'][idx] = dPD
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
    results['sigma_z'][idx] = sigma_z
    results['uw'][idx] = uw
    results['fR_bathurst'][idx] = fR_bathurst
    results['fR_colebrook_white'][idx] = fR_colebrook_white
    
    results['hw'][idx] = hw
    results['SCs'][idx] = SCs
    results['Qout'][idx] = Qout
    results['Qin_compensated'][idx] = Qin_compensated
    results['Vadd_E'][idx] = Vadd_E
    results['Vadd_C'][idx] = Vadd_C
    results['Vadd_TM'][idx] = Vadd_TM
    
    
    
'''Plot'''
# pmm.plot_geom(results, time, z, set_xlim_moulin=False,
#               ax2_varname=['Mr_major','Mr_minor'],
#               ax3_varname=['dC_major','dC_minor'],
#               ax4_varname=['dE_major','dE_minor'], 
#               ax5_varname='dTM',              
#               ax6_varname='dPD' )
#pmm.plot_geom(results,ax2_varname=['dC_major','dC_minor'],ax3_varname='Pw_z',ax4_varname='Pi_z',ax5_varname='sigma_z'  )
#pmm.plot_2Darray(results,results['dTM'])


'''Plot water level and cross-section area timeseries'''
# pmm.plot_1Darray_timeserie(results, time, 'hw')
# pmm.plot_1Darray_timeserie(results, time, 'SCs')


'''Plots volumes'''
# plt.figure()
# #plt.plot(time,Qin,label='Qin')
# plt.plot(time,results['Vadd_C'],label='Vadd_C')
# plt.plot(time,results['Vadd_E'],label='Vadd_E')
# plt.plot(time,results['Vadd_TM'],label='Vadd_TM')
# plt.legend()


# plt.figure()
# plt.subplot(2,1,1)
# Vadd_C_percent = results['Vadd_C'].flatten()/Qin *100
# Vadd_E_percent = results['Vadd_E'].flatten()/Qin *100
# plt.plot(time, Vadd_C_percent, label='Vadd_C')
# plt.plot(time, Vadd_E_percent, label='Vadd_E')
# plt.title('Percent Qin volume change')
# plt.legend()
# plt.grid(True)
# plt.subplot(2,1,2)
# plt.plot(time,Qin,label='Qin')
# plt.plot(time,results['Qin_compensated'].flatten(),label='Qin compensated')
# plt.legend()
# plt.grid(True)





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


# #%%
# x_lim=10

# fig,ax = plt.subplots()
# ax.axhspan(0, hw, facecolor ='lightblue', alpha = 1,zorder=1)
# ax.axhspan(-100, 0, facecolor ='peru', alpha = 1,zorder=1)

# ax.plot(Mx_upstream,z,color='black') #plot major axis on the left
# ax.plot(Mx_downstream,z,color='black')  #plot minor axis on the right
# ax.set_xlim([-x_lim,x_lim])


# ax.fill_betweenx(z,-x_lim,Mx_upstream,color='aliceblue',zorder=2)
# ax.fill_betweenx(z,Mx_downstream,x_lim,color='aliceblue',zorder=2)
