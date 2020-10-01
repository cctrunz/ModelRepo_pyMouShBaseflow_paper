import function_moulin_model as fmm
import matplotlib.pyplot as plt


Mr_top=2
Mr_bottom=2

H = 1000
hw = 800
[z, nz, dz] = fmm.generate_grid_z(H) 
[Mx_upstream, Mx_downstream, Mr_major, Mr_minor]= fmm.initiate_moulin_wall(
    z,type='linear',Mr_top=Mr_top,Mr_bottom=Mr_bottom)


sigma_x = 0e3#-50e3 #(Units??) compressive
sigma_y = 50e3#-50e3 #(Units??) compressive
tau_xy = -50e3#100e3 #(Units??) shear opening


wet = fmm.locate_water(hw,z) 
Pw_z = fmm.calculate_water_pressure_at_depth(hw,z,wet) #water pressure at each depth
Pi_z = fmm.calculate_ice_pressure_at_depth(H,z) #ice pressure at each depth
sigma_z = fmm.calculate_sigma_z(Pw_z, Pi_z)


[dE_major,dE_minor] = fmm.calculate_elastic_deformation(
    Mr_major, Mr_minor, sigma_z, sigma_x, sigma_y, tau_xy)

plt.figure()
plt.plot(dE_major,z)

