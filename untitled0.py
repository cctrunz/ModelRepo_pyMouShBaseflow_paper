import numpy as np
import matplotlib.pyplot as plt
import function_moulin_model as fmm


dt = 300 #(s)
hw1 = [200,300,400,500,600,700,800,900,1000]
hw2 = [500,600,700,800,900,950,1000]
H = 1000
z = np.arange(H)
E = 5
Mr = np.linspace(4,1,len(z))
T_mean = np.linspace(273.15,273.15,len(z))
Pi_z=fmm.calculate_ice_pressure_at_depth(H, z)
iceflow_param_glen = fmm.calculate_iceflow_law_parameter(T_mean,Pi_z)

''' Exploration impact of creep'''
    
plt.figure()
plt.plot(Pi_z,z,color='red',label='Ice pressure')
plt.plot(0,0,color='blue',label='Water pressure')
plt.plot([0,1e7],[900,900],color='black',label='Overburden pressure')
#plt.plot(Pw_z,z,color='red')
#plt.plot(fmm.calculate_water_pressure_at_depth(hw,z),z,color='black')
for i in hw1:   
    wet = fmm.locate_water(i,z) 
    Pw_z=fmm.calculate_water_pressure_at_depth(i,z,wet)
    plt.plot(Pw_z,z,color='blue')
plt.legend()
plt.ylabel('z')
plt.xlabel('Pressure')

plt.figure()
plt.plot([0,0],[0,1000],color='black')
for i in hw2:   
    wet = fmm.locate_water(i,z) 
    Pw_z=fmm.calculate_water_pressure_at_depth(i,z,wet)
    plt.plot((Pw_z-Pi_z)/1000,z,linestyle='--')
plt.legend()
plt.ylabel('z')
plt.xlabel('Pressure')



#Mr = 2*np.ones(len(z))

plt.figure(figsize=[8,5])
for i in hw2: 
    wet = fmm.locate_water(i,z) 
    Pw_s=fmm.calculate_water_pressure_at_depth(i,z,wet)   
    creep = fmm.calculate_creep_moulin(Mr,Mr,dt,iceflow_param_glen,-Pi_z,Pw_z,E)
    plt.plot(creep[1],z)
plt.ylabel('z')
plt.xlabel('Creep (m)')

plt.figure(figsize=[3,5])
plt.plot(Mr,z,-Mr,z,color='black')
plt.xlim([-5,5])
#plt.title('Moulin radius')