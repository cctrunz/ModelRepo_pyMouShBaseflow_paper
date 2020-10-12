import function_moulin_model as fmm
import matplotlib.pyplot as plt
import pandas as pd
#import numpy as np

temperature_profil_litterature = pd.read_excel('FieldData/Temp_Profil_Litterature.xlsx',sheet_name=None) #import all the sheet at the same time
print(temperature_profil_litterature.keys()) #list the temperature profiles options
temperature_profile = temperature_profil_litterature['IkenB'].temperature #[273.15,273.15]
temperature_profile = temperature_profile.iloc[::-1]#reverse data to have the 0 m elevation z match the 0 row
# plt.figure()
# plt.plot(temperature_profile)

H = 1000
hw = 600
Mr_top=1
Mr_bottom=1

dz = 1
z = fmm.generate_grid_z(H,dz) 
Mr = fmm.initiate_moulin_radius(z,type='linear',Mr_top=Mr_top,Mr_bottom=Mr_bottom)

lim_x = [-2,2]

#%% dE

sigma_x = 0e3#-50e3 #(Units??) compressive
sigma_y = 0e3#-50e3 #(Units??) compressive
tau_xy = 0e3#100e3 #(Units??) shear opening

wet = fmm.locate_water(hw,z) 
Pw_z = fmm.calculate_water_pressure_at_depth(hw,z,wet) #water pressure at each depth
Pi_z = fmm.calculate_ice_pressure_at_depth(H,z) #ice pressure at each depth
sigma_z = fmm.calculate_sigma_z(Pw_z, Pi_z)

dE = fmm.calculate_elastic_deformation(Mr, sigma_z, sigma_x, sigma_y, tau_xy)
plt.figure()
plt.plot([0,0],[0,H],'--',color='black')
plt.plot(dE*1000,z,label='dE')

plt.legend()
plt.xlabel('(mm)')
plt.ylabel('z')
plt.xlim(lim_x)
plt.savefig('Figures/d_save/dE.pdf')


#%% dC

dt = 300
E = 1
T_ice = fmm.interpolate_T_profile(z,temperature_profile=[273.15,273.15])#temperature_profile)
wet = fmm.locate_water(hw,z) 
Pw_z = fmm.calculate_water_pressure_at_depth(hw,z,wet) #water pressure at each depth
Pi_z = fmm.calculate_ice_pressure_at_depth(H,z) #ice pressure at each depth
sigma_z = fmm.calculate_sigma_z(Pw_z, Pi_z)
iceflow_param_glen = fmm.calculate_iceflow_law_parameter(T_ice,Pi_z)

dC = fmm.calculate_creep_moulin(Mr,dt,iceflow_param_glen,sigma_z,E)
plt.figure()
plt.plot([0,0],[0,H],'--',color='black')
plt.plot(dC*1000*60*60/300,z,label='dC')
plt.legend()
plt.xlabel('(mm/h)')
plt.ylabel('z')
plt.xlim(lim_x)
plt.savefig('Figures/d_save/dC.pdf')

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
T_ice = fmm.interpolate_T_profile(z,temperature_profile=[273.15,273.15])#temperature_profile)

dTM = fmm.calculate_melt_below_head(dL,head_loss_dz_TM, dt, Qout, Mpr, wet,include_ice_temperature=True,T_ice=T_ice,Tmw=Tmw)
plt.figure()
plt.plot([0,0],[0,H],'--',color='black')
plt.plot(dTM*1000*60*60/300,z,label='dTM')
plt.legend()
plt.xlabel('(mm/h)')
plt.ylabel('z')
plt.xlim(lim_x)
plt.savefig('Figures/d_save/dM.pdf')



plt.figure()
plt.plot(T_ice,z,label='$T_{ice}$')
plt.legend()
plt.xlabel('(degree Farenheit)')
plt.ylabel('z')
plt.savefig('Figures/d_save/Tice.pdf')


# #%% dOC

# Qin = 2
# dt = 300
# friction_factor_OC = 0.5

# wet = fmm.locate_water(hw,z) 
# dL = fmm.calculate_dL(Mr, dz)
# [Mcs, Mpr, Mdh, Mrh] = fmm.calculate_moulin_geometry(Mr, Mr)
# uw_OC = fmm.calculate_water_velocity(Qin, Mcs)#,wet)
# head_loss_dz_OC = fmm.calculate_moulin_head_loss(uw_OC, friction_factor_OC, dL, Mdh)
# T_ice = fmm.interpolate_T_profile(z,temperature_profile=temperature_profile)#[273.15,273.15])

# dOC = fmm.calculate_melt_above_head_OC(Mr,Mpr,dL,dt,head_loss_dz_OC,Qin,wet,include_ice_temperature=True,T_ice=T_ice)
# plt.figure()
# plt.plot([0,0],[0,H],'--',color='black')
# plt.plot(dOC*1000*60*60/300,z,label='dOC')
# plt.legend()
# plt.xlabel('(mm/h)')
# plt.ylabel('z')
# plt.xlim(lim_x)
# #plt.savefig('Figures/d_save/dOC.pdf')




# #%% dF

# t=300
# dt = 300
# wet = fmm.locate_water(hw,z) 
# T_ice = fmm.interpolate_T_profile(z,temperature_profile=[260,260])

# dF = fmm.calculate_refreezing(T_ice,t,dt,wet)
# dF_test = fmm.calculate_refreezing_test(T_ice,t,dt,wet)

# plt.figure()
# plt.plot([0,0],[0,H],'--',color='black')
# plt.plot(dF*1000*60*60/300,z,label='dF')
# plt.legend()
# plt.xlabel('(mm/h)')
# plt.ylabel('z')
# plt.xlim(lim_x)
# #plt.savefig('Figures/d_save/dF.pdf')

# plt.figure()
# plt.plot([0,0],[0,H],'--',color='black')
# plt.plot(dF_test*1000*60*60/300,z,label='dF test')
# plt.legend()
# plt.xlabel('(mm/h)')
# plt.ylabel('z')
# plt.xlim(lim_x)
# #plt.savefig('Figures/d_save/dFtest.pdf')

# #%%

