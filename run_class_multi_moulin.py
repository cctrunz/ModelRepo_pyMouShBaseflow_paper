from PyMouSh import MoulinShape, TimeStamps, Qin_constant, Qin_sinusoidal, Qin_real, calculate_h_S_schoof
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#%%
secinday = 24*3600
ZERO_KELVIN = 273.15


jeme_basin = pd.read_csv('Field_Data/surface_melt_jeme.csv')
jeme_basin = jeme_basin.dropna()
Qin_data = jeme_basin.Qm3s.to_numpy() +0.1
Qtime_data = jeme_basin.SOY.to_numpy()

jeme_moulin = pd.read_csv('Field_Data/head_jeme.csv')
jeme_moulin = jeme_moulin.dropna()
h_real = jeme_moulin.head_bed.to_numpy()
t_real = jeme_moulin.soy.to_numpy()

time_start = Qtime_data[0]
time_end = Qtime_data[0] + 10*secinday#Qtime_data[-1]#time_start + 5*secinday #
timestep = 300 #seconds
time = TimeStamps(time_start,time_end,timestep)

meltwater_input_moulin2 = Qin_real(time, Qin_data, Qtime_data)
meltwater_input_moulin1 = meltwater_input_moulin2 *10
# Qin_mean = 1
# dQ = 0.1 
# meltwater_input = Qin_sinusoidal(time,Qin_mean, dQ)


#paramters
                                
dz = 1
z_elevations = None
moulin_radii = 0.5
temperature_profile=np.array([ZERO_KELVIN, ZERO_KELVIN])
ice_thickness = 500
initial_head = 500
initial_subglacial_area = 1      
regional_surface_slope = 0
channel_length = 500
creep_enhancement_factor = 5
sigma = (0, 0)  # -50e3 #(Units??) compressive #50e3#-50e3
tau_xy = 0  # -50e3#100e3 #(Units??) shear opening
friction_factor_OC = 1
friction_factor_TM = 0.5
friction_factor_SUB = 0.1
fraction_pd_melting = 0.2



moulin1 = MoulinShape(                      
                    dz = dz,
                    z_elevations = z_elevations,
                    moulin_radii = moulin_radii,                  
                    temperature_profile = temperature_profile,                   
                    ice_thickness = ice_thickness,
                    initial_head = initial_head,
                    initial_subglacial_area = initial_subglacial_area,
                            
                    regional_surface_slope = regional_surface_slope,
                    channel_length = channel_length,
                    creep_enhancement_factor = creep_enhancement_factor,
                    
                    sigma = sigma,  # -50e3 #(Units??) compressive #50e3#-50e3
                    tau_xy = tau_xy,  # -50e3#100e3 #(Units??) shear opening
                    friction_factor_OC = friction_factor_OC,
                    friction_factor_TM = friction_factor_TM,
                    friction_factor_SUB = friction_factor_SUB,
                    fraction_pd_melting = fraction_pd_melting)


moulin2 = MoulinShape(                      
                    dz = dz,
                    z_elevations = z_elevations,
                    moulin_radii = moulin_radii,                  
                    temperature_profile = temperature_profile,                   
                    ice_thickness = ice_thickness,
                    initial_head = initial_head,
                    initial_subglacial_area = initial_subglacial_area,
                            
                    regional_surface_slope = regional_surface_slope,
                    channel_length = channel_length,
                    creep_enhancement_factor = creep_enhancement_factor,
                    
                    sigma = sigma,  # -50e3 #(Units??) compressive #50e3#-50e3
                    tau_xy = tau_xy,  # -50e3#100e3 #(Units??) shear opening
                    friction_factor_OC = friction_factor_OC,
                    friction_factor_TM = friction_factor_TM,
                    friction_factor_SUB = friction_factor_SUB,
                    fraction_pd_melting = fraction_pd_melting)



head_L_moulin1 = None 

for idx,t in enumerate(time) :
    #main subglacial channel
    moulin1.run1step(time,
                    timestep,
                    meltwater_input_moulin1[idx],
                    subglacial_baseflow = 0, 
                    head_L = head_L_moulin1 )
    #jeme    
    moulin2.run1step(time,
                    timestep,
                    meltwater_input_moulin2[idx],
                    subglacial_baseflow = 3, 
                    head_L = moulin1.head )    
    
    
    
    
    
    #print(moulin.head)
   
#moulin.listdict[0].keys()    
#moulin.dict.keys()








#%%
# fig, ax = plt.subplots()
# idx = -2
# moulin1.plot_head(ax)
# moulin2.plot_head(ax)
plt.figure()
plt.plot(moulin1.dict['head'],label='moulin 1')
plt.plot(moulin2.dict['head'],label='moulin 2')
plt.legend()

#%%

#%%

# plt.figure()
# plt.plot(time,moulin.dict['meltwater_output_subglacial'] )

# plt.figure()
# plt.plot(time,moulin.dict['head'] )

# plt.figure()
# for idx in np.arange(0,len(time),100):
#     plt.plot(moulin.listdict[idx]['moulin_wall_position_upstream'],moulin.z)
#     plt.plot(moulin.listdict[idx]['moulin_wall_position_downstream'],moulin.z)
#     plt.xlim([-2,2])
#     plt.pause(0.1)
#%%

# def plot_AGU(moulin, 
#               time,
#               idx,
#               t_real,
#               h_real,
#               spine_head_min=200
#               ):

#     time = time/24/3600
#     z = moulin.z
#     ice_thickness = moulin.ice_thickness    
#     head = moulin.dict['head']
#     subglacial_area = moulin.dict['subglacial_cross_section_area']   
#     meltwater_input_moulin = moulin.dict['meltwater_input_moulin']
#     #melwater_input_compensated = moulin.dict['melwater_input_compensated_moulin']
#     meltwater_output_subglacial = moulin.dict['meltwater_output_subglacial'] 
#     Mx_upstream = moulin.listdict[idx]['moulin_wall_position_upstream']
#     Mx_downstream = moulin.listdict[idx]['moulin_wall_position_downstream']        
      
#     fig = plt.figure(figsize=(25,8))
#     grid = plt.GridSpec(4,4)#, wspace=-0.7)
#     ax1 = fig.add_subplot(grid[0:4, 0])
#     ax2 = fig.add_subplot(grid[0:4, 1:4])#, sharey=ax1)  #hw
#     ax3 = fig.add_subplot(grid[2, 1:4])#ax2.twinx() #SCs
#     ax4 = fig.add_subplot(grid[3, 1:4], sharex=ax2)    #Qin-Qout
#     ax5 = ax4.twinx()#fig1.add_subplot(grid[2, 1:4], sharex=ax2)    #Qin-Qout
    
#     # # dlim = 10
#     # # dbound = 10
#     # # dtick = [-10,5,0,5,10]
    
#     time_plot=time[0:idx+1]
#     min_time = int(min(time))
#     max_time = int(max(time))
    
#     ax1.clear() #clear figure content -- especially important for ploting in the loop
#     ax2.clear()
#     ax3.clear()
#     ax4.clear()
#     ax5.clear() 
    
#     # #MOULIN SHAPE
    
#     ax1.plot(Mx_upstream,z,color='black') #plot major axis on the left
#     ax1.plot(Mx_downstream,z,color='black')  #plot minor axis on the right
    
#     ax2.plot(t_real/3600/24,h_real,'-',color='black')      
#     ax2.plot(time_plot,head[0:idx+1],'-',color='blue')   
#     ax2.legend(['measured','simulated'],loc="upper left")
     
#     ax3.plot(time_plot,subglacial_area[0:idx+1],'-',color='red')       
    
#     ax4.plot(time,meltwater_input_moulin,'--',color='blue')#,label='Qin')    
#     ax5.plot(time_plot,meltwater_output_subglacial[0:idx+1],color='red')#,'-',color='red',label='Qout')
    
#     ax1.axhspan(0, head[idx], facecolor ='lightblue', alpha = 1,zorder=1)
#     ax1.axhspan(-140, 0, facecolor ='peru', alpha = 1,zorder=1)
#     ax1.fill_betweenx(z,-10,Mx_upstream,color='aliceblue',zorder=2)
#     ax1.fill_betweenx(z,Mx_downstream,10,color='aliceblue',zorder=2)     
    
#     '''Moulin'''
#     ax1.set_ylabel('z(m)')    
#     ax1.set_xlabel('(m)')
#     ax1.set_xlim([-5,5]) 
#     ax1.set_ylim([-50,ice_thickness]) 
#     ax1.spines['top'].set_visible(False)
#     ax1.spines['right'].set_visible(False)
#     ax1.spines['bottom'].set_position(('zero')) 
#     ax1.spines['left'].set_bounds(0,ice_thickness)
#     ax1.spines['bottom'].set_bounds(-5,5)
#     ax1.set_xticks([-4,-3,-2,-1,0,1,2,3,4]) 
#     ax1.set_xticklabels([-4,-3,-2,-1,0,1,2,3,4]) 
#     ax1.set_yticks(np.arange(0,ice_thickness+1,100)) 
#     ax1.set_yticklabels(np.arange(0,ice_thickness+1,100))
    
#     '''Head'''
#     ax2.set_ylabel('Head',color='blue')
#     ax2.set_xlim([min_time,max_time])
#     ax2.set_ylim([-50,ice_thickness]) 
#     ax2.yaxis.tick_left()
#     ax2.yaxis.set_label_position("left")
#     ax2.tick_params(axis='y', labelcolor='blue')
#     ax2.spines['top'].set_visible(False)
#     ax2.spines['bottom'].set_visible(False)
#     ax2.spines['right'].set_visible(False)
#     ax2.spines['left'].set_color('blue')
#     ax2.axes.xaxis.set_visible(False)
#     ax2.spines['left'].set_bounds(spine_head_min,ice_thickness)
#     ax2.set_yticks(np.arange(spine_head_min,ice_thickness+1,100)) 
#     ax2.set_yticklabels(np.arange(spine_head_min,ice_thickness+1,100))
    
    
#     '''Subglacial channel'''
#     ax3.set_ylabel('Subglacial_area',color='red')
#     ax3.set_xlim([min_time,max_time])
#     ax3.set_ylim([min(subglacial_area),max(subglacial_area)])
#     ax3.yaxis.tick_right()
#     ax3.yaxis.set_label_position("right")
#     ax3.tick_params(axis='y', labelcolor='red')
#     ax3.spines['top'].set_visible(False)
#     ax3.spines['bottom'].set_visible(False)
#     ax3.spines['left'].set_visible(False)
#     ax3.spines['right'].set_color('red')
#     ax3.axes.xaxis.set_visible(False)
    
    
#     '''Qin'''
#     ax4.set_ylabel('Qin',color='blue')
#     ax4.set_xlim([min_time,max_time])
#     ax4.set_ylim([min(meltwater_input_moulin),max(meltwater_input_moulin)])
#     ax4.spines['top'].set_visible(False)
#     ax4.spines['right'].set_visible(False)
#     ax4.spines['left'].set_color('blue')
#     ax4.tick_params(axis='y', labelcolor='blue')
#     ax4.set_xticks(np.round(np.linspace(min_time,20,max_time)))
#     ax4.set_xticklabels(np.round(np.linspace(min_time,20,max_time)))
    
#     '''Qout'''
#     ax5.set_ylabel('Qout',color='red')
#     ax5.set_xlim([min_time,max_time])
#     ax5.set_ylim([min(meltwater_output_subglacial),max(meltwater_output_subglacial)])
#     ax5.yaxis.tick_right()
#     ax5.yaxis.set_label_position("right")
#     ax5.tick_params(axis='y', labelcolor='red')
#     ax5.spines['top'].set_visible(False)
#     ax5.spines['left'].set_visible(False)
#     ax5.spines['bottom'].set_visible(False)
#     ax5.spines['right'].set_color('red')
    
    
    
#     #make backgroud transparent    
#     ax1.patch.set_alpha(0)
#     ax2.patch.set_alpha(0)
#     ax3.patch.set_alpha(0)
#     ax4.patch.set_alpha(0)
#     ax5.patch.set_alpha(0)
     
#%%
# idx = 3000
# plot_AGU(moulin,time,idx,t_real,h_real,spine_head_min=200)   
                   

# for idx in np.arange(0,len(time/timestep),10):
#     plot_AGU(moulin,time,idx,t_real,h_real,spine_head_min=200)   
#     plt.savefig('')
    
   
    
   
    