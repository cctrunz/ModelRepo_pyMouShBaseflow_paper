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
#import plot_codes.plot_moulin_model as pmm
import pandas as pd
#import plot_codes.comprehensive_plot

# import plot_codes.plot_pretty_moulin
# import plot_codes.plot_deltas_all
sim_number = 4
param_filename = 'Saved_Sim/Parameters_Sim_%d'%sim_number
constants_filename = 'Saved_Sim/Constants_Sim_%d'%sim_number
results_filename = 'Saved_Sim/Results_Sim_%d'%sim_number


'''Glacier parameters'''
#R0 = 2 #(m) Initial moulin radius
H = 500 #(m) Ice thickness
regional_surface_slope =0# fmm.calculate_alpha(H)#0.01#alpha in matlab %regional surface slope (unitless), for use in Glen's Flow Law
L = 15000#fmm.calculate_L(H) #L = 10000 #(m) Subglacial channel length 
E = 0.1 #Enhancement factor for the ice creep.
#Assign elastic deformation parameters
sigma_x = 0#-50e3 #(Units??) compressive
sigma_y = 0#50e3#-50e3 #(Units??) compressive
tau_xy = 0#-50e3#100e3 #(Units??) shear opening

#Turbulent melting parameters
friction_factor_OC = 0.1
friction_factor_TM = 0.5
fraction_pd_melting = 0.01

'''Model parameters'''
# tmax_in_day = 10 #(days) Maximum time to run
# dt = 300 #(s) timestep
# time = fmm.generate_time(dt,tmax_in_day)
dz = 1 #(m) 

z = fmm.generate_grid_z(H,dz) #default is 1m spacing # should this mimic the x grid??

'''Qin from the field'''
# baseflow = 0.10 #m3/s
# Qin_filename ='lowc17_Q_jeme_rolling'#'lowc18_Q_pira_rolling'#high17_Q_radi_rolling,high18_Q_radi_rolling,lowc17_Q_jeme_rolling
# Qin_array = pd.read_csv('Melt_field_data/smooth_melt_data/%s.csv'%Qin_filename,parse_dates=True, index_col='Date') + baseflow
# tmax_in_day = 120#Qin_array.Seconds[len(Qin_array)-1] /3600/24 #10 #(days) Maximum time to run
# dt = 300 #(s) timestep
# time = fmm.generate_time(dt,tmax_in_day)
# Qin = fmm.set_Qin(time,type='field_data', Qin_array=Qin_array.melt_rate, time_array=Qin_array.Seconds)#*120

'''Qin fake'''
tmax_in_day = 50
dt = 300 #(s) timestep
time = fmm.generate_time(dt,tmax_in_day)
Qin_mean = 2 #m3/s
period = 24*3600
dQ = 0.1
Qin = fmm.set_Qin(time,type='sinusoidal_celia',Qin_mean=Qin_mean,dQ=dQ,period=period)
Qin_filename = 'mean1_dQ0.1'
baseflow = '–'


'''Moulin parameters'''
Mr_top=0.2
Mr_bottom=0.2
Mr_minimum = 1e-9 #(m)

Mr_major = fmm.initiate_moulin_radius(z,type='linear',Mr_top=Mr_top,Mr_bottom=Mr_bottom)
# z_custom = [0,299,300,H]
# Mr_custom = [0.2,0.2,0.2,0.2]#[4,4,0.2,0.2]
# Mr_major = np.interp(z,z_custom,Mr_custom)
Mr_minor = Mr_major
[Mx_upstream, Mx_downstream]= fmm.initiate_moulin_wall_position(Mr_major,Mr_minor) #for position of the wall relative to coordinate system

'''Temperature profile'''
#import some temperature profiles from the litterature											
temperature_profil_litterature = pd.read_excel('FieldData/Temp_Profil_Litterature.xlsx',sheet_name=None) #import all the sheet at the same time
print(temperature_profil_litterature.keys()) #list the temperature profiles options
temperature_profile = temperature_profil_litterature['Lüthi15_FOXX1'].temperature #[273.15,273.15]
temperature_profile = temperature_profile.iloc[::-1]#reverse data to have the 0 m elevation z match the 0 row
#temperature_profile = [273.15,273.15]
T_ice = fmm.interpolate_T_profile(z,temperature_profile=temperature_profile)

'''Calculate parameters (out of the loop)'''
#---calculate total ice pressure
Pi_H = fmm.calculate_ice_pressure_at_depth(H,0) #total ice pressure (ice pressure at the bottom)
Pi_z = fmm.calculate_ice_pressure_at_depth(H,z) #ice pressure at each depth
iceflow_param_glen = fmm.calculate_iceflow_law_parameter(T_ice,Pi_z) #(units?) called A in matlab's code

secinday=24*3600
mts_to_cmh = 100*60*60/dt #m per timestep to mm/h : change units

'''Initial values'''
hw = H #(m)Initial water level
SCs = 1 #(m) Initial subglacial channel croohhhss-section area
#initialize
dGlen = 0
dGlen_cumulative = 0
Vadd_C = 0
Vadd_E = 0
Vadd_TM = 0
#Initiate results dictionnary
results = fmm.initiate_results_dictionnary(time,z)

idx_plot = 0
idx_save = 0
i_save = 0


param = {'H':H,'L':L,'E':E,'alpha':regional_surface_slope,\
             'sigma_x':sigma_x,'sigma_y':sigma_y,'tau_xy':tau_xy,\
             'fr_factor_OC':friction_factor_OC,'fr_factor_TM':friction_factor_TM,'fraction_pd_melting':fraction_pd_melting,\
             'z':z,'time':time,'Mr_initial':Mr_major,
             'Qin':Qin,'T_ice':T_ice, 'Qin_filename':Qin_filename, 'baseflow':baseflow, \
             'hw_initial':hw, 'SCs_initial':SCs,\
             'dE':'on','dC':'on','dM':'on','dPD':'on','dOC':'off'
             }
fmm.pickle_dictionnary(param,param_filename)

constants = fmm.save_constants()
fmm.pickle_dictionnary(constants,constants_filename)



for idx, t in enumerate(time):
      
    '''Calculate moulin geometry, relative to water level'''    
    [Mcs, Mpr, Mdh, Mrh] = fmm.calculate_moulin_geometry(Mr_major, Mr_minor)
    
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
    uw_TM = fmm.calculate_water_velocity(Qout, Mcs)
    uw_OC = fmm.calculate_water_velocity(Qin[idx], Mcs)
    #friction_factor = fmm.calculate_relative_friction_factor(Mdh,Mrh,relative_roughness,type='unknown')
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
    
    #Potential drop
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
    dr_major = fmm.calculate_dradius(dE=dE_major, dC=dC_major, dTM=dTM, dPD=dPD)
    dr_minor = fmm.calculate_dradius(dE=dE_minor, dC=dC_minor, dTM=dTM, dPD=dPD)
    Mr_major = fmm.update_moulin_radius( Mr_major,dr_major )
    Mr_minor = fmm.update_moulin_radius( Mr_minor,dr_minor )
    # if any(Mr_major)<=0:
    #     Mr_major.all[Mr_major<=0]=0.05
    # if any(Mr_minor)<=0:
    #     Mr_minor.all[Mr_minor<=0]=0.05
    [Mx_upstream, Mx_downstream] = fmm.update_moulin_wall_position(Mx_upstream, Mx_downstream, dr_major,dr_minor, dGlen, dGlen_cumulative)
    

    '''Save values'''
    results['Mx_upstream'][idx] = Mx_upstream
    results['Mx_downstream'][idx] = Mx_downstream
    results['Mr_major'][idx] = Mr_major
    results['Mr_minor'][idx] = Mr_minor
    results['dC_major'][idx] = dC_major
    results['dC_minor'][idx] = dC_minor
    results['dTM'][idx] = dTM
    #results['dOC'][idx] = dOC
    results['dGlen'][idx] =  dGlen
    results['dGlen_cumulative'][idx] =  dGlen_cumulative
    results['dE_major'][idx] = dE_major
    results['dE_minor'][idx] = dE_minor
    results['dr_major'][idx] = dr_major
    results['dr_minor'][idx] = dr_minor
    results['Mcs'][idx] =  Mcs
    results['Mpr'][idx] = Mpr
    results['Mdh'][idx] = Mdh
    results['Mrh'][idx] = Mrh
    results['Pi_z'][idx] = Pi_z
    results['Pw_z'][idx] = Pw_z
    results['wet'][idx] = wet
    results['Tmw'][idx] = Tmw
    results['sigma_z'][idx] = sigma_z
    results['uw_TM'][idx] = uw_TM
    results['uw_OC'][idx] = uw_OC
    results['friction_factor_TM'][idx] = friction_factor_TM
    results['friction_factor_OC'][idx] = friction_factor_OC

    
    results['hw'][idx] = hw
    results['SCs'][idx] = SCs
    results['Qout'][idx] = Qout
    results['Qin_compensated'][idx] = Qin_compensated
    results['Vadd_E'][idx] = Vadd_E
    results['Vadd_C'][idx] = Vadd_C
    results['Vadd_TM'][idx] = Vadd_TM

fmm.pickle_dictionnary(results,results_filename)
    
#     if idx_plot == idx:
#          idx_plot = idx_plot+20
#          idx_save = idx_save+1
#     #     #plot_codes.plot_deltas_all.live_plot(dTM,dOC,dOC,dC_major,dC_minor,dE_major,dE_minor,dr_major,dr_minor,dGlen,z,mts_to_cmh)
#     #     plot_codes.plot_pretty_moulin.live_plot(hw,Mx_upstream,Mx_downstream,z) 
#          plt.pause(0.001)
# #plot_codes.plot_pretty_moulin.live_plot(hw,Mx_upstream,Mx_downstream,z)   
#     #plot_codes.plot_pretty_moulin.live_plot(hw,Mx_upstream,Mx_downstream,z) 
#     #plt.savefig('Pira_moulin_first_attempt/Pira_moulin_first_attempt_%s'%idx_save)
    
#          plot_codes.comprehensive_plot_new.live_plot(Mx_upstream,Mx_downstream,dTM,dC_major,dC_minor,dE_major,dE_minor,dr_major,dr_minor,dGlen,t,hw,SCs,z,Qin,Qout,idx,wet,mts_to_cmh,time,H)
        


# colors = [plt.cm.rainbow(i) for i in np.linspace(0, 1, len(time))] 
# plt.figure()
# for i in np.arange(len(time)):
#     #plt.plot(results['Mr_major'][i],results['z']) 
#     plt.plot(results['dOC'][i],z,color=colors[i]) 
# plt.figure()
# plt.plot(time,Qin)
# plt.figure()
# plt.plot(time,results['hw'])

#%%

'''Import field data for comparison'''
jeme = pd.read_csv('Head_field_data/JEME17_MOU_1.csv',parse_dates=True, index_col='TIMESTAMP')
radi = pd.read_csv('Head_field_data/RADI17_MOU_1.csv',parse_dates=True, index_col='TIMESTAMP')
pira = pd.read_csv('Head_field_data/PIRA18_HYD_1.csv',parse_dates=True, index_col='TIMESTAMP')

'''create elapsed time array'''
elapsed = jeme.index-jeme.index[0] 
jeme['Seconds']=elapsed.total_seconds() +1.56555e+06
elapsed = radi.index-radi.index[0] 
radi['Seconds']=elapsed.total_seconds() +83700
elapsed = pira.index-pira.index[0] 
pira['Seconds']=elapsed.total_seconds()+ 163800




'''Comprehensive plot for movies'''
# import plot_codes.comprehensive_plot_new_with_results
# plot_interval = 20
# plot_codes.comprehensive_plot_new_with_results.live_plot(results,dt,z,Qin,time,H,plot_interval)

'''Single head-moulin plots'''
import plot_codes.moulin_and_head_plot_with_results
plot_codes.moulin_and_head_plot_with_results.live_plot(results,dt,z,Qin,time,H,jeme.Seconds,jeme.water_level_above_bed,idx=idx)
plt.savefig('Saved_Sim/Plot_Sim_%d'%sim_number)

#%%
'''Loop head-moulin plots'''
# import plot_codes.moulin_and_head_plot_with_results

# idx_plot = 0
# idx_save = 0
# i_save = 0
# for idx,t in enumerate(time):
#     if idx_plot == idx:
#           idx_plot = idx_plot+5
#           idx_save = idx_save+1
#           plot_codes.moulin_and_head_plot_with_results.live_plot(results,dt,z,Qin,time,H,jeme.Seconds,jeme.water_level_above_bed,idx=idx)
#           plt.pause(0.001)
#           #plt.savefig('Movies/Figure_movie_%s'%idx_save)





