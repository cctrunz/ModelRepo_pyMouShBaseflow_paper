''' 
Moulin Physical Model from LC Andrews and K Poinar
Translated by Celia Trunz

Component provenance:
    - Subglacial channel: Schoof2010
    - 
'''

import numpy as np
import function_moulin_model as fmm


'''Default parameter'''
R0 = 2 #(m) Initial moulin radius
H = 500 #(m) Ice thickness
tmax_in_day = 20 #(days) Maximum time to run
dt = 300 #(s) timestep
r_minor_initial = 2
r_major_initial = 2
xmax    = 30 # 80 #(m) how far away from moulin to use as infinity
hw = H/2 #Initial water level


#Assign elastic deformation parameters
stress={} #initiate dictionnary for elastic deformation parameters
stress['sigma_x'] = -50e3 #(Units??) compressive
stress['sigma_y'] = -50e3 #(Units??) compressive
stress['tau_xy'] = 100e3 #(Units??) shear opening

#
delta={} #initiate dictionnary for moulin shape variation

''' set model grid'''
z = fmm.generate_grid_z(H) #default is 1m spacing
#M ={} #create and initialize the M, containing the moulin characterisitcs #use this to have all the moulin data in one dictionnary

Pi = fmm.ice_pressure_at_depth(H,0) #total ice pressure

''' set moulin geom '''
[Mr_minor, Mr_major] = fmm.initiate_moulin_radius(z,r_minor_initial,r_major_initial) #
[x, nx, dx] = fmm.generate_grid_x(dt, xmax)

''' set duration of the model run'''
[time,tmax_in_second] = fmm.generate_time(dt,tmax_in_day)

'''set ice temperature'''
[T_far,T] = fmm.set_ice_temperature(x,z)


for t in time:
   # wet = locate_water(hw,z)
    
    Pi_z = fmm.ice_pressure_at_depth(H,z) #ice pressure at each depth
    Pw_z = fmm.water_pressure_at_depth(hw,z)
    
    wet = fmm.locate_water(hw,z) 
    stress['cryo'] = -Pi_z # Ice hydrostatic stress (INWARD: Negative)
    stress['hydro'] = Pw_z # Water hydrostatic stress (OUTWARD: Positive)'
    stress['hydro'][np.invert(wet)] = 0 # There is no stress from the water above the water level
    

    delta['creep_moulin_minor'] = fmm.creep_moulin(Mr_minor,dt,T,Pi_z,stress)
   # print(delta_creep_moulin_minor)