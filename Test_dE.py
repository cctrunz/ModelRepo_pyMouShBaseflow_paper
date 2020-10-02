import function_moulin_model as fmm
import matplotlib.pyplot as plt


Mr_top=2
Mr_bottom=2

H = 1000
hw = 800
dz = 1
z = fmm.generate_grid_z(H,dz) 
Mr = fmm.initiate_moulin_radius(
    z,type='linear',Mr_top=Mr_top,Mr_bottom=Mr_bottom)



sigma_x = 0e3#-50e3 #(Units??) compressive
sigma_y = 50e3#-50e3 #(Units??) compressive
tau_xy = -50e3#100e3 #(Units??) shear opening


wet = fmm.locate_water(hw,z) 
Pw_z = fmm.calculate_water_pressure_at_depth(hw,z,wet) #water pressure at each depth
Pi_z = fmm.calculate_ice_pressure_at_depth(H,z) #ice pressure at each depth
sigma_z = fmm.calculate_sigma_z(Pw_z, Pi_z)

dt = 300 #(s) timestep
xmax    = 30
[x, dx] = fmm.generate_grid_x(dt, xmax)
[T_far,T] = fmm.set_ice_temperature(x,z,ice_temperature=[273.15,273.15])


dE = fmm.calculate_elastic_deformation(
    Mr, sigma_z, sigma_x, sigma_y, tau_xy)

plt.figure()
plt.plot(dE,z)

