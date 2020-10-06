import function_moulin_model as fmm
import matplotlib.pyplot as plt
import numpy as np

Mr_top=1
Mr_bottom=1

H = 1000
hw = 800
S = 2
L = 10000
E = 5
dz = 1
dt = 300
Qin = 2

friction_factor_OC = 0.5
friction_factor_TM = 0.1

z = fmm.generate_grid_z(H,dz) 
Mr = fmm.initiate_moulin_radius(
    z,type='linear',Mr_top=Mr_top,Mr_bottom=Mr_bottom)

[Mcs, Mpr, Mdh, Mrh] = fmm.calculate_moulin_geometry(Mr, Mr)

sigma_x = 0e3#-50e3 #(Units??) compressive
sigma_y = 50e3#-50e3 #(Units??) compressive
tau_xy = -50e3#100e3 #(Units??) shear opening


wet = fmm.locate_water(hw,z) 
Pw_z = fmm.calculate_water_pressure_at_depth(hw,z,wet) #water pressure at each depth
Pi_z = fmm.calculate_ice_pressure_at_depth(H,z) #ice pressure at each depth
sigma_z = fmm.calculate_sigma_z(Pw_z, Pi_z)
Tmw = fmm.calculate_pressure_melting_temperature(Pw_z) 

T_ice = fmm.interpolate_T_profile(z,temperature_profile=[273.15,273.15])
iceflow_param_glen = fmm.calculate_iceflow_law_parameter(T_ice,Pi_z)


Qout = Qin#fmm.calculate_Qout(S,hw,L)



uw_TM = fmm.calculate_water_velocity(Qout, Mcs)#,wet)
uw_OC = fmm.calculate_water_velocity(Qin, Mcs)#,wet)
dL = fmm.calculate_dL(Mr, dz)
head_loss_dz_TM = fmm.calculate_moulin_head_loss(uw_TM, friction_factor_TM, dL, Mdh)
head_loss_dz_OC = fmm.calculate_moulin_head_loss(uw_OC, friction_factor_OC, dL, Mdh)

dE = fmm.calculate_elastic_deformation(Mr, sigma_z, sigma_x, sigma_y, tau_xy)
dC = fmm.calculate_creep_moulin(Mr,dt,iceflow_param_glen,sigma_z,E)
dTM = fmm.calculate_melt_below_head(dL,head_loss_dz_TM, dt, Qout, Mpr, wet,include_ice_temperature=True,T_ice=T_ice,Tmw=Tmw)
dOC = fmm.calculate_melt_above_head_OC(Mr,Mpr,dL,dt,head_loss_dz_OC,Qin,wet,include_ice_temperature=False,T_ice=T_ice)

# T0 = 273.15 #Kelvin = 0°C
# rhoi = 910 #kg/m3; Ice density
# rhow = 1000 #kg/m3; Water density
# g = 9.8 #m/s2; Gravity
# ki = 2.1 #J/mKs
# cp = 2115 #J/kgK
# Lf = 335000 #J/kg; Latent heat of fusion
# cw = 4210

# remove_neg = np.zeros(len(dL))
# remove_neg[dL>=0]=1    
# AA = rhow * g * Qin * (head_loss_dz_OC /dL)
# BB = Mpr * rhoi * (cw *(T0-T_ice)+Lf)

# dOC = AA/BB
# dOC[wet]=0
# dOC_dt=dOC*dt*remove_neg




plt.figure()
# plt.plot(dE,z)
# plt.plot(dC,z)
plt.plot(dTM,z)
plt.plot(dOC,z)


#%%

# plt.figure()
# for radius in np.linspace(0.1,5,5):
#     Mr_top=radius
#     Mr_bottom=radius
#     Mr = fmm.initiate_moulin_radius(z,type='linear',Mr_top=Mr_top,Mr_bottom=Mr_bottom)
#     [Mcs, Mpr, Mdh, Mrh] = fmm.calculate_moulin_geometry(Mr, Mr)

#     dL = fmm.calculate_dL(Mr, dz)
#     head_loss_dz_TM = fmm.calculate_moulin_head_loss(uw, friction_factor_TM, dL, Mdh)
#     dTM = fmm.calculate_melt_below_head(dL,head_loss_dz_TM, dt, Qout, Mpr, wet,include_ice_temperature=True,T_ice=T_ice,Tmw=Tmw)
#     plt.semilogx(dTM,z,label='r=%s' %radius)
# plt.legend()
# plt.xlabel('Turbulent melt')
# plt.ylabel('z')

# #%%

# T0 = 273.15 #Kelvin = 0°C
# rhoi = 910 #kg/m3; Ice density
# rhow = 1000 #kg/m3; Water density
# g = 9.8 #m/s2; Gravity
# cp = 2115 #J/kgK
# Lf = 335000 #J/kg; Latent heat of fusion
# cw = 4210 



# # remove_neg = np.zeros(len(dL))
# # remove_neg[dL>=0]=1    

# dOC_dt = (rhow * g * Qin * (head_loss_dz_OC /dL)) / (Mpr * rhoi * (cw *(T0-T_ice)+Lf))
    
# dOC_dt[wet]=0
# # dOC_dt=dOC_dt*remove_neg
# plt.figure()
# plt.plot(dOC_dt,z)
    