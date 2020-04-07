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



'''Default parameter'''
#R0 = 2 #(m) Initial moulin radius
H = 1000 #(m) Ice thickness
L = 10000 #(m) Subglacial channel length
tmax_in_day = 2 #(days) Maximum time to run
dt = 300 #(s) timestep
r_minor_initial = 1 #(m)
r_major_initial = 1 #(m)
xmax    = 30 # 80 #(m) how far away from moulin to use as infinity
hw = H/2 #(m)Initial water level
S = 0.5 #(m)
Qin = 2 #m3/s


#Assign elastic deformation parameters
stress={} #initiate dictionnary for elastic deformation parameters
stress['sigma_x'] = -50e3 #(Units??) compressive
stress['sigma_y'] = -50e3 #(Units??) compressive
stress['tau_xy'] = 100e3 #(Units??) shear opening

#Turbulent melting parameters
relative_roughness = 0.2 #increasing this value increases the amount of melting due to turbulence. (comment from Matlab code)
relative_roughness_OC = 1e-9 #1e-12;  % This one modifies the melt from open channel flow.
#
#delta={} #initiate dictionnary for moulin shape variation
#stress={}


''' set model grid'''
[z, nz, dz] = fmm.generate_grid_z(H) #default is 1m spacing # should this mimic the x grid??
#M ={} #create and initialize the M, containing the moulin characterisitcs #use this to have all the moulin data in one dictionnary

Pi = fmm.ice_pressure_at_depth(H,0) #total ice pressure

''' set moulin geom '''
[Mr_minor, Mr_major] = fmm.initiate_moulin_radius(z,r_minor_initial,r_major_initial) #
[x, nx, dx] = fmm.generate_grid_x(dt, xmax) #generate grid around the moulin for temperature discretization

''' set duration of the model run'''
[time,tmax_in_second] = fmm.generate_time(dt,tmax_in_day)

'''set ice temperature'''
[T_far,T] = fmm.set_ice_temperature(x,z)

plt.figure()
plt.plot(Mr_minor,z)

for idx, t in enumerate(time):
    
    '''Calculate or import meltwater input'''
    Qin = fmm.calculate_Qin(t,type='constant',Qin_mean=3)
    
    '''Calculate moulin geometry'''
    [Ms, Mp, Mdh, Mrh] = fmm.calculate_moulin_geometry(Mr_minor, Mr_major)
    
    '''Calculate water level''' 
    sol = solve_ivp(fmm.function_subglacial_schoof,
                    [0, dt], #initial time and end time. We solve for one timestep.
                    [hw,S], #initial head and channel cross-section area. Uses values in previous timestep.
                    args=(Ms,z,Pi,L,Qin), #args() only works with scipy>=1.4. if below, then error message: missing 5 arguments
                    method = 'LSODA' #solver method
                    # atol = 1e-6, #tolerance. can make it faster
                    # rtol = 1e-3,
                    #max_step = 10 #change if resolution is not enough
                    )
    hw = sol.y[0][-1]  #(m) moulin water level
    S = sol.y[1][-1] #(m) Channel cross-section
    Qout = fmm.calculate_Qout(S,hw,L)


    '''Calculate parameters'''
    Pi_z = fmm.ice_pressure_at_depth(H,z) #ice pressure at each depth
    Pw_z = fmm.water_pressure_at_depth(hw,z) #water pressure at each depth
    wet = fmm.locate_water(hw,z) 
    stress_cryo = -Pi_z # Ice hydrostatic stress (INWARD: Negative)
    stress_hydro = Pw_z # Water hydrostatic stress (OUTWARD: Positive)'
    stress_hydro[np.invert(wet)] = 0 # There is no stress from the water above the water level
    
    
    '''Calculate moulin changes for each component'''
    #Creep Deformation
    delta_creep_moulin_minor = fmm.calculate_creep_moulin(Mr_minor,dt,T,Pi_z,stress_cryo,stress_hydro)
    delta_creep_moulin_major = fmm.calculate_creep_moulin(Mr_major,dt,T,Pi_z,stress_cryo,stress_hydro)
    
    #Turbulent melting
    
    
    
    #Refreezing
    
    
    #Open channel melting
    
    #Elastic deformation
    
    #Asymmetric deformation due to Glen's Flow Law
    
    
    
    '''Calculate new moulin radii'''
    
    Mr_minor = Mr_minor + delta_creep_moulin_minor
    Mr_major = Mr_major + delta_creep_moulin_major
    
    #
    
    #print(delta['creep_moulin_minor'])
    
    
    #plt.plot(delta_creep_moulin_major,z)
    #plt.plot(t,hw,'.',color='red')
    #plt.plot(t,Qout,'.',color='black')
    plt.plot(Mr_minor,z)