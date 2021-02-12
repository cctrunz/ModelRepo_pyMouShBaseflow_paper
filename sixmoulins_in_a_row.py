from PyMouSh import MoulinShape, TimeStamps, Qin_constant, Qin_sinusoidal, Qin_real, calculate_h_S_schoof
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle

#%%
secinday = 24*3600
ZERO_KELVIN = 273.15

temperature_profile=np.array([ZERO_KELVIN, ZERO_KELVIN])
   
# friction_factor_OC = 1
# friction_factor_TM = 0.5
# friction_factor_SUB = 0.1
# fraction_pd_melting = 0.2

jeme_moulin = pd.read_csv('Field_Data/head_jeme.csv')
jeme_moulin = jeme_moulin.dropna()
h_real = jeme_moulin.head_bed.to_numpy()
t_real = jeme_moulin.soy.to_numpy()

jeme_basin = pd.read_csv('Field_Data/surface_melt_jeme.csv')
jeme_basin = jeme_basin.dropna()
Qin_data = jeme_basin.Qm3s.to_numpy() + 0.1
Qtime_data = jeme_basin.SOY.to_numpy() + 3*3600 #add routing delay -- Jess found 2h not 4. investigate why

time_start = Qtime_data[int(2*secinday/900)]  
time_end = time_start + 50*secinday
timestep = 300 #seconds
time = TimeStamps(time_start,time_end,timestep)

meltwater_input1 = Qin_real(time, Qin_data, Qtime_data)
meltwater_input2 = Qin_real(time + 1*3600, Qin_data, Qtime_data)
meltwater_input3 = Qin_real(time + 2*3600, Qin_data, Qtime_data)
meltwater_input4 = Qin_real(time + 3*3600, Qin_data, Qtime_data)
meltwater_input5 = Qin_real(time + 4*3600, Qin_data, Qtime_data)
meltwater_input6 = Qin_real(time + 5*3600, Qin_data, Qtime_data)
meltwater_input7 = Qin_real(time + 6*3600, Qin_data, Qtime_data)

distance_to_margin = [5000,9000,13000,17000,21000,25000]

H0= 500./np.sqrt(25000.) #assures ice thickness of 500 m at 25 km from edge
ice_thickness = H0*np.sqrt(distance_to_margin)
initial_head = ice_thickness

moulin1 = MoulinShape(ice_thickness = ice_thickness[0],
                      initial_head = initial_head[0],
                      initial_subglacial_area = np.pi*0.25**2/2, 
                      channel_length = 5e3)

moulin2 = MoulinShape(ice_thickness = ice_thickness[1],
                      initial_head = initial_head[1],
                      initial_subglacial_area = np.pi*0.25**2/2, 
                      channel_length = 5e3)

moulin3 = MoulinShape(ice_thickness = ice_thickness[2],
                      initial_head = initial_head[2],
                      initial_subglacial_area = np.pi*0.25**2/2, 
                      channel_length = 5e3)

moulin4 = MoulinShape(ice_thickness = ice_thickness[3],
                      initial_head = initial_head[3],
                      initial_subglacial_area = np.pi*0.25**2/2, 
                      channel_length = 5e3)

moulin5 = MoulinShape(ice_thickness = ice_thickness[4],
                      initial_head = initial_head[4],
                      initial_subglacial_area = np.pi*0.5**2/2, 
                      channel_length = 5e3)

moulin6 = MoulinShape(ice_thickness = ice_thickness[5],
                      initial_head = initial_head[5],
                      initial_subglacial_area = np.pi*1**2/2, 
                      channel_length = 5e3)




for idx,t in enumerate(time) :
    #margin moulin
    moulin1.run1step(time,
                    timestep,
                    meltwater_input1[idx],
                    subglacial_baseflow = 0, 
                    head_L = None )
       
    moulin2.run1step(time,
                    timestep,
                    meltwater_input2[idx],
                    subglacial_baseflow = 0,  
                    head_L = moulin1.head )  
    
    moulin3.run1step(time,
                    timestep,
                    meltwater_input3[idx],
                    subglacial_baseflow = 0,  
                    head_L = moulin2.head )    
        
    moulin4.run1step(time,
                    timestep,
                    meltwater_input4[idx],
                    subglacial_baseflow = 0,  
                    head_L = moulin3.head )    
            
    moulin5.run1step(time,
                    timestep,
                    meltwater_input5[idx],
                    subglacial_baseflow = 0,  
                    head_L = moulin4.head )    
    
    moulin6.run1step(time,
                    timestep,
                    meltwater_input6[idx],
                    subglacial_baseflow = 0,  
                    head_L = moulin5.head )    
    
    
picklefile = open('moulin1', 'wb')
pickle.dump(moulin1, picklefile)
picklefile.close()

picklefile = open('moulin2', 'wb')
pickle.dump(moulin2, picklefile)
picklefile.close()

picklefile = open('moulin3', 'wb')
pickle.dump(moulin3, picklefile)
picklefile.close()

picklefile = open('moulin4', 'wb')
pickle.dump(moulin4, picklefile)
picklefile.close()

picklefile = open('moulin5', 'wb')
pickle.dump(moulin5, picklefile)
picklefile.close()

picklefile = open('moulin6', 'wb')
pickle.dump(moulin6, picklefile)
picklefile.close()


#%% plot all the head timeserie for all the moulins

time_year = time/secinday
fig, axs = plt.subplots(4, sharex=False, figsize=(7,7))
fig.tight_layout() 

x_profile = np.arange(25000)
ice_thickness_profile =H0*np.sqrt(x_profile)
axs[0].plot(x_profile,ice_thickness_profile, color = 'black')
axs[0].plot([distance_to_margin[0],distance_to_margin[0]],[ice_thickness[0],0], linewidth=5, label='moulin 1')
axs[0].plot([distance_to_margin[1],distance_to_margin[1]],[ice_thickness[1],0], linewidth=5, label='moulin 2')
axs[0].plot([distance_to_margin[2],distance_to_margin[2]],[ice_thickness[2],0], linewidth=5, label='moulin 3')
axs[0].plot([distance_to_margin[3],distance_to_margin[3]],[ice_thickness[3],0], linewidth=5, label='moulin 4')
axs[0].plot([distance_to_margin[4],distance_to_margin[4]],[ice_thickness[4],0], linewidth=5, label='moulin 5')
axs[0].plot([distance_to_margin[5],distance_to_margin[5]],[ice_thickness[5],0], linewidth=5, label='moulin 6')
axs[0].set_ylabel('Ice thickness (m)')
axs[0].set_xlabel('Distance from margin (m)')
axs[0].invert_xaxis()
axs[0].legend(loc='upper right')

axs[1].plot(time_year,moulin1.dict['meltwater_input_moulin'])
axs[1].plot(time_year,moulin2.dict['meltwater_input_moulin'])
axs[1].plot(time_year,moulin3.dict['meltwater_input_moulin'])
axs[1].plot(time_year,moulin4.dict['meltwater_input_moulin'])
axs[1].plot(time_year,moulin5.dict['meltwater_input_moulin'])
axs[1].plot(time_year,moulin6.dict['meltwater_input_moulin'])
axs[1].set_ylabel('Qin ($m^3/s$)')
axs[1].set_xlim([180,200])

axs[2].plot(t_real/secinday,h_real, color='black')
axs[2].plot(time_year, moulin1.dict['head'],label='moulin 1')
axs[2].plot(time_year, moulin2.dict['head'],label='moulin 2')
axs[2].plot(time_year, moulin3.dict['head'],label='moulin 3')
axs[2].plot(time_year, moulin4.dict['head'],label='moulin 4')
axs[2].plot(time_year, moulin5.dict['head'],label='moulin 5')
axs[2].plot(time_year, moulin6.dict['head'],label='moulin 6')
axs[2].set_ylabel('Head (m)')
axs[2].set_xlim([180,200])
axs[2].set_ylim([0,500])

axs[3].plot(time_year, moulin1.dict['subglacial_radius'],label='moulin 1')
axs[3].plot(time_year, moulin2.dict['subglacial_radius'],label='moulin 2')
axs[3].plot(time_year, moulin3.dict['subglacial_radius'],label='moulin 3')
axs[3].plot(time_year, moulin4.dict['subglacial_radius'],label='moulin 4')
axs[3].plot(time_year, moulin5.dict['subglacial_radius'],label='moulin 5')
axs[3].plot(time_year, moulin6.dict['subglacial_radius'],label='moulin 6')
axs[3].set_ylabel('sub radius (m)')
axs[3].set_xlabel('day of year')
axs[3].set_xlim([180,200])


#%% plot the moulin position in the ice sheet

plt.figure()



