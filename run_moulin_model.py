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

'''Activate or deactivate False and True statements'''

include_ice_temperature = False #%true means that the change in the ice temperature is included in...
#%the calculated change in moulin radius. If false, it makes the implicit
#%assumption that the ice temperature and water temperature are both at the pressure melting temperature.

#do the same for Bathurst....




'''Default parameter'''
#R0 = 2 #(m) Initial moulin radius
H = 1000 #(m) Ice thickness
L = 10000 #(m) Subglacial channel length
tmax_in_day = 10 #(days) Maximum time to run
dt = 300 #(s) timestep
Mr_minor_initial = 2 #(m)
Mr_major_initial = 2 #(m)
Mr_minimum = 1e-9 #(m)
xmax    = 30 # 80 #(m) how far away from moulin to use as infinity
hw = H #(m)Initial water level
Ss = 1.5 #(m) Initial subglacial channel cross-section area
Qin = 3 #m3/s
E = 3 #Enhancement factor for the ice creep.

#Assign elastic deformation parameters
stress={} #initiate dictionnary for elastic deformation parameters
stress['sigma_x'] = -50e3 #(Units??) compressive
stress['sigma_y'] = -50e3 #(Units??) compressive
stress['tau_xy'] = 100e3 #(Units??) shear opening

#Turbulent melting parameters
relative_roughness = 1 #increasing this value increases the amount of melting due to turbulence. (comment from Matlab code)
relative_roughness_OC = 1e-9 #1e-12;  % This one modifies the melt from open channel flow.

'''Initialize arrays'''
#delta={} #initiate dictionnary for moulin shape variation
#stress={}
#M ={} #create and initialize the M, containing the moulin characterisitcs #use this to have all the moulin data in one dictionnary

'''set and initialize parameters and dimensions'''
#---set model grid
[z, nz, dz] = fmm.generate_grid_z(H) #default is 1m spacing # should this mimic the x grid??
#---set moulin geom 
[Mr_major, Mr_minor]= fmm.initiate_moulin_radius(z,Mr_major_initial,Mr_minor_initial) #for moulin dimension
#[Mx_major, Mx_minor]= fmm.initiate_moulin_wall_position(Mr_major, Mr_minor) #for position of the wall relative to coordinate system
[x, nx, dx] = fmm.generate_grid_x(dt, xmax) #generate grid around the moulin for temperature discretization
#---set duration of the model run'
[time,tmax_in_second] = fmm.generate_time(dt,tmax_in_day)
#---set ice temperature
[T_far,T] = fmm.set_ice_temperature(x,z)
#specific for turbulence
if include_ice_temperature:
    Ti = T_far
else:
    Ti = np.nan

'''Calculate parameters'''
#---calculate total ice pressure
Pi = fmm.ice_pressure_at_depth(H,0) 


#Initiate results dictionnary
results={}
results['Mr_major']= np.zeros([len(time),len(z)])
results['Mr_minor'] = np.zeros([len(time),len(z)])
results['dC_major'] = np.zeros([len(time),len(z)]) 
results['dC_minor'] = np.zeros([len(time),len(z)]) 
results['dTM_major'] = np.zeros([len(time),len(z)]) 
results['dTM_minor'] = np.zeros([len(time),len(z)]) 
results['Ms'] = np.zeros([len(time),len(z)]) 
results['Mp'] = np.zeros([len(time),len(z)]) 
results['Mdh'] = np.zeros([len(time),len(z)]) 
results['Mrh'] = np.zeros([len(time),len(z)]) 
results['Pi_z'] = np.zeros([len(time),len(z)]) 
results['Pw_z'] = np.zeros([len(time),len(z)]) 
results['wet'] = np.zeros([len(time),len(z)]) 
results['Pw'] = np.zeros([len(time),len(z)]) 
results['Tmw'] = np.zeros([len(time),len(z)]) 
results['stress_cryo'] = np.zeros([len(time),len(z)]) 
results['stress_hydro'] = np.zeros([len(time),len(z)]) 
results['uw'] = np.zeros([len(time),len(z)]) 
results['fR_bathurst'] = np.zeros([len(time),len(z)])
results['fR_colebrook_white'] = np.zeros([len(time),len(z)]) 

results['hw'] = np.zeros([len(time),1])
results['Ss'] = np.zeros([len(time),1])
results['Qin'] = np.zeros([len(time),1])
results['Qout'] = np.zeros([len(time),1])
results['time'] = time
results['z'] = z
#save also initial radius and parameters??

#makes an artificial cone
# Mr_major=Mr_major*np.linspace(1,2,len(z))
# Mr_minor=Mr_minor*np.linspace(1,2,len(z))







for idx, t in enumerate(time):
    
    '''Calculate or import meltwater input'''
    #sinusoidal_celia
    #constant
    Qin = fmm.calculate_Qin(t,type='sinusoidal_celia',Qin_mean=3)
    
    '''Calculate moulin geometry'''
    [Ms, Mp, Mdh, Mrh] = fmm.calculate_moulin_geometry(Mr_minor, Mr_major)
    
    '''Calculate water level''' 
    sol = solve_ivp(fmm.function_subglacial_schoof,
                    [0, dt], #initial time and end time. We solve for one timestep.
                    [hw,Ss], #initial head and channel cross-section area. Uses values in previous timestep.
                    args=(Ms,z,Pi,L,Qin), #args() only works with scipy>=1.4. if below, then error message: missing 5 arguments
                    method = 'LSODA' #solver method
                    # atol = 1e-6, #tolerance. can make it faster
                    # rtol = 1e-3,
                    #max_step = 10 #change if resolution is not enough
                    )
    hw = sol.y[0][-1]  #(m) moulin water level
    Ss = sol.y[1][-1] #(m) Channel cross-section
    Qout = fmm.calculate_Qout(Ss,hw,L)


    '''Calculate parameters'''
    Pi_z = fmm.ice_pressure_at_depth(H,z) #ice pressure at each depth
    Pw_z = fmm.water_pressure_at_depth(hw,z) #water pressure at each depth
    wet = fmm.locate_water(hw,z) 
    Tmw=fmm.calculate_pressure_melting_temperature(Pw_z)   
    stress_cryo = -Pi_z # Ice hydrostatic stress (INWARD: Negative)
    stress_hydro = Pw_z # Water hydrostatic stress (OUTWARD: Positive)'
    stress_hydro[np.invert(wet)] = 0 # There is no stress from the water above the water level
    uw = fmm.calculate_water_velocity(Qout, Ms, wet)
    fR_bathurst = fmm.calculate_bathurst_friction_factor(Mrh, relative_roughness)
    fR_colebrook_white = fmm.calculate_colebrook_white_friction_factor(Mdh, relative_roughness)
    
    
    
    '''Calculate moulin changes for each component'''
    #Creep Deformation

    dC_major = fmm.calculate_creep_moulin(Mr_major,dt,T,Pi_z,stress_cryo,stress_hydro,E)
    dC_minor = fmm.calculate_creep_moulin(Mr_minor,dt,T,Pi_z,stress_cryo,stress_hydro,E)
    
    #Turbulent melting
    [dTM_major,Vadd_turb] = fmm.calculate_turbulent_melting_moulin(
        Mr_major, 0.01, uw, Tmw, Ti, z, dz, dt, Qout, Mp, Mdh, include_ice_temperature)
    [dTM_minor,Vadd_turb] = fmm.calculate_turbulent_melting_moulin(
        Mr_minor, 0.01, uw, Tmw, Ti, z, dz, dt, Qout, Mp, Mdh, include_ice_temperature)
    
    
    #Refreezing
    
    
    #Open channel melting
    
    #Elastic deformation
    
    #Asymmetric deformation due to Glen's Flow Law
    
    
    
    '''Calculate new moulin radii'''
    
    [Mr_major, Mr_minor] = fmm.calculate_new_moulin_radius(
        Mr_major, Mr_minor, 
        dC=[dC_major, dC_minor],
        dTM=[dTM_major,dTM_minor])

    
    
    #Mr_minor = Mr_minor + delta_creep_moulin_minor
    #Mr_major = Mr_major + delta_creep_moulin_major
    
    #
    
    #print(delta['creep_moulin_minor'])
    
    '''Save values'''
    results['Mr_major'][idx] = Mr_major
    results['Mr_minor'][idx] = Mr_minor
    results['dC_major'][idx] = dC_major
    results['dC_minor'][idx] = dC_minor
    results['dTM_major'][idx] = dTM_major
    results['dTM_minor'][idx] = dTM_minor
    results['Ms'][idx] =  Ms
    results['Mp'][idx] = Mp
    results['Mdh'][idx] = Mdh
    results['Mrh'][idx] = Mrh
    results['Pi_z'][idx] = Pi_z
    results['Pw_z'][idx] = Pw_z
    results['wet'][idx] = wet
    results['Tmw'][idx] = Tmw
    results['stress_cryo'][idx] = stress_cryo
    results['stress_hydro'][idx] = stress_hydro
    results['uw'][idx] = uw
    results['fR_bathurst'][idx] = fR_bathurst
    results['fR_colebrook_white'][idx] = fR_colebrook_white
    
    results['hw'][idx] = hw
    results['Ss'][idx] = Ss
    results['Qin'][idx] = Qin
    results['Qout'][idx] = Qout
    
    
    
    #plt.plot(delta_creep_moulin_major,z)
    #plt.plot(t,hw,'.',color='red')
    #plt.plot(t,Qout,'.',color='black')
    #plt.plot(Mr_minor,z)


#%%
plt.figure()
colors = [plt.cm.rainbow(i) for i in np.linspace(0, 1, len(results['time']))] 
for i in np.arange(0,len(results['time']),100):
#for i in np.arange(len(results['time'])):
    plt.plot(-results['Mr_major'][i],results['z'],color=colors[i])
    plt.plot(results['Mr_major'][i],results['z'],color=colors[i])
    #plt.plot(results['dTM_major'][i],results['z'])
    #plt.xlim([0,1e-9])
    #plt.plot(results['dC_major'][i],results['z'],color='blue')

plt.figure()    
for i in np.arange(0,len(results['time']),100):
#for i in np.arange(len(results['time'])):
    plt.plot(results['dC_major'][i],results['z'],color=colors[i])

plt.figure()
for i in np.arange(0,len(results['time']),100):
#for i in np.arange(len(results['time'])):
    plt.plot(results['dTM_major'][i],results['z'],color=colors[i])    
    plt.title('Melting')
#%%
#plt.imshow(np.rot90(results['dTM_major']))#,origin='lower'
#plt.colorbar()
#plt.clim(1e-9,1.2e-9)

