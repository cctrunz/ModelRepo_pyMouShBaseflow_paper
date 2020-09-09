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

'''Activate or deactivate False and True statements'''









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
E = 5 #Enhancement factor for the ice creep.
regional_surface_slope = 0.02#alpha in matlab %regional surface slope (unitless), for use in Glen's Flow Law
#Q_type = 'constant' #
Q_type = 'sinusoidal_celia'
#import some temperature profiles from the litterature											
temperature_profil_litterature = pd.read_excel('FieldData/Temp_Profil_Litterature.xlsx',sheet_name=None) #import all the sheet at the same time
print(temperature_profil_litterature.keys()) #list the temperature profiles options
Temperature_profile = temperature_profil_litterature['Lüthi15_FOXX1'].temperature #[273.15,273.15]
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
relative_roughness = 1 #increasing this value increases the amount of melting due to turbulence. (comment from Matlab code)
relative_roughness_OC = 1e-9 #1e-12;  % This one modifies the melt from open channel flow.

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
#Save parameter in dictionnary


#makes an artificial cone
# Mr_major=Mr_major*np.linspace(1,2,len(z))
# Mr_minor=Mr_minor*np.linspace(1,2,len(z))


fig = plt.figure(figsize=(20,5))
#fig.patch.set_alpha(1)
# ax1 = fig.add_subplot(131)
# ax2 = fig.add_subplot(194,sharey=ax1)
# ax3 = fig.add_subplot(195,sharey=ax1)
# ax4 = fig.add_subplot(196,sharey=ax1)
# ax5 = fig.add_subplot(197,sharey=ax1)
# ax6 = fig.add_subplot(198,sharey=ax1)
# plt.subplots_adjust(wspace=-0.2)

grid = plt.GridSpec(1, 16, wspace=-0.7)
ax1 = fig.add_subplot(grid[0, 1:4])
ax2 = fig.add_subplot(grid[0, 5:8], sharey=ax1)
ax3 = fig.add_subplot(grid[0, 7:10], sharey=ax1)
ax4 = fig.add_subplot(grid[0, 9:12], sharey=ax1)
ax5 = fig.add_subplot(grid[0, 11:14], sharey=ax1)
ax6 = fig.add_subplot(grid[0, 13:16], sharey=ax1)




idx_plot = 0
#plt.ion()
for idx, t in enumerate(time):
      
    '''Calculate moulin geometry, relative to water level'''
    
    [Mcs, Mpr, Mdh, Mrh,Diameter] = fmm.calculate_moulin_geometry(
        Mx_upstream, Mx_downstream, Mr_major, Mr_minor)
    
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


    '''Calculate parameters (in loop)'''


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
    dTM = fmm.calculate_melt_below_head(
        Mx_upstream, Mx_downstream, fR_bathurst, uw, Tmw, Ti, z, dz, dt, Qout, Mpr, Mdh, include_ice_temperature,wet)

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
    
    
    '''Calculate new moulin radii'''
 
    
    [Mx_upstream, Mx_downstream, Mr_major, Mr_minor] = fmm.calculate_new_moulin_wall_position( Mx_upstream, Mx_downstream, dGlen_cumulative,
                                                                    dC = [dC_major, dC_minor],
                                                                    dTM = dTM,
                                                                    dGlen = dGlen,
                                                                    dE = [dE_major,dE_minor],
                                                                    dPD = dPD
                                                                    )
    dr_major = dTM + dPD + dC_major + dE_major
    dr_minor = dTM + dPD + dC_minor + dE_minor
    
    if idx_plot == idx:
        #
        idx_plot = idx_plot+20
        # pmm.plot_pretty_moulin(Mx_upstream,Mx_downstream,hw,z,x_lim=10)
        # plt.pause(0.1)
        
        ax1.clear() #clear figure content -- especially important for ploting in the loop
        ax1.plot(Mx_upstream,z,color='black') #plot major axis on the left
        ax1.plot(Mx_downstream,z,color='black')  #plot minor axis on the right
                  
        ax2.clear()
        ax2.plot(dTM[wet]*mts_to_cmh,z[wet],color='black') #plot major axis on the left
        ax2.plot(dPD[~wet]*mts_to_cmh,z[~wet],color='black') #plot major axis on the left
        
        ax3.clear()
        ax3.plot(dC_major*mts_to_cmh,z,color='grey') #plot major axis on the left
        ax3.plot(dC_minor*mts_to_cmh,z,color='black')  #plot minor axis on the right

        ax4.clear()
        ax4.plot(dE_major*mts_to_cmh,z,color='grey') #plot major axis on the left
        ax4.plot(dE_minor*mts_to_cmh,z,color='black')  #plot minor axis on the right  
        
        ax5.clear()
        ax5.plot(dr_major*mts_to_cmh,z,color='grey') #plot major axis on the left
        ax5.plot(dr_minor*mts_to_cmh,z,color='black')  #plot minor axis on the right  
        
        ax6.clear()
        ax6.plot(dGlen*mts_to_cmh,z,color='black') #plot major axis on the left
        
        
        ax2.set_title('Turbulent Melting')
        ax3.set_title('Creep Deformation')
        ax4.set_title('Elastic Deformation')
        ax5.set_title('Change in Radius')
        ax6.set_title('Ice Motion')
        
        ax1.set_ylabel('z(m)')
        
        ax1.set_xlabel('(m)')
        ax2.set_xlabel('(cm/h)')
        ax3.set_xlabel('(cm/h)')
        ax4.set_xlabel('(cm/h)')
        ax5.set_xlabel('(cm/h)')
        ax6.set_xlabel('(cm/h)')

        ax1.set_xlim([-10,10]) 
        ax2.set_xlim([-5,5])
        ax3.set_xlim([-5,5])
        ax4.set_xlim([-5,5])
        ax6.set_xlim([-5,5])
        ax5.set_xlim([-5,5])
        
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
    
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        
        ax3.spines['top'].set_visible(False)
        ax3.spines['right'].set_visible(False)
        ax3.spines['left'].set_visible(False)

        ax4.spines['top'].set_visible(False)
        ax4.spines['right'].set_visible(False)
        ax4.spines['left'].set_visible(False)
        
        ax5.spines['top'].set_visible(False)
        ax5.spines['right'].set_visible(False)
        ax5.spines['left'].set_visible(False)
        
        ax6.spines['top'].set_visible(False)
        ax6.spines['right'].set_visible(False)
        ax6.spines['left'].set_visible(False)
        
        ax2.axes.yaxis.set_visible(False)
        ax3.axes.yaxis.set_visible(False)
        ax4.axes.yaxis.set_visible(False)
        ax5.axes.yaxis.set_visible(False)
        ax6.axes.yaxis.set_visible(False)
        
        ax1.spines['bottom'].set_position(('zero'))
        ax2.spines['bottom'].set_position(('zero'))
        ax3.spines['bottom'].set_position(('zero'))
        ax4.spines['bottom'].set_position(('zero'))
        ax5.spines['bottom'].set_position(('zero'))
        ax6.spines['bottom'].set_position(('zero'))
        
        ax1.spines['left'].set_bounds(0,H)
        ax1.spines['bottom'].set_bounds(-10,10)
        ax2.spines['bottom'].set_bounds(-1,1)
        ax3.spines['bottom'].set_bounds(-1,1)
        ax4.spines['bottom'].set_bounds(-1,1)
        ax5.spines['bottom'].set_bounds(-1,1)
        ax6.spines['bottom'].set_bounds(-1,1)
    
        ax1.set_xticks([-8,-6,-4,-2,0,2,4,6,8]) 
        ax2.set_xticks([-1,0,1]) 
        ax3.set_xticks([-1,0,1])
        ax4.set_xticks([-1,0,1])
        ax5.set_xticks([-1,0,2])
        ax6.set_xticks([-1,0,1])
        ax1.set_xticklabels([8,6,4,2,0,2,4,6,8]) 
        ax2.set_xticklabels([-1,0,1])         
        ax3.set_xticklabels([-1,0,1])         
        ax4.set_xticklabels([-1,0,1])         
        ax5.set_xticklabels([-1,0,1])         
        ax6.set_xticklabels([-1,0,1])
        
        #ax3.patch.set_facecolor('red')
        ax2.patch.set_alpha(0)
        ax2.patch.set_alpha(0)
        ax3.patch.set_alpha(0)
        ax4.patch.set_alpha(0)
        ax5.patch.set_alpha(0)
        ax6.patch.set_alpha(0)
        
        ax1.axhspan(0, hw, facecolor ='lightblue', alpha = 1,zorder=1)
        ax1.axhspan(-140, 0, facecolor ='peru', alpha = 1,zorder=1)
        ax1.fill_betweenx(z,-10,Mx_upstream,color='aliceblue',zorder=2)
        ax1.fill_betweenx(z,Mx_downstream,10,color='aliceblue',zorder=2)        
        ax2.fill_betweenx(z[wet],0,dTM[wet]*mts_to_cmh,color='orangered')
        ax2.fill_betweenx(z[~wet],0,dPD[~wet]*mts_to_cmh,color='orangered')
        ax3.fill_betweenx(z,dC_minor*mts_to_cmh,0,where=dC_minor>0,facecolor=('orangered'))
        ax3.fill_betweenx(z,dC_minor*mts_to_cmh,0,where=dC_minor<0,facecolor=('lightgreen'))
        ax4.fill_betweenx(z,dE_minor*mts_to_cmh,0,where=dE_minor>0,facecolor=('orangered'))
        ax4.fill_betweenx(z,dE_minor*mts_to_cmh,0,where=dE_minor<0,facecolor=('lightgreen'))
        ax5.fill_betweenx(z,dr_minor*mts_to_cmh,0,where=dr_minor>0,facecolor=('orangered'))
        ax5.fill_betweenx(z,dr_minor*mts_to_cmh,0,where=dr_minor<0,facecolor=('lightgreen'))
        ax6.fill_betweenx(z,0,dGlen*mts_to_cmh,color='grey')
        plt.pause(0.1)


    

    #plt.pause(0.5)
    
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
