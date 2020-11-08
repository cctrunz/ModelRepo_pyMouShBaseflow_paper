from PyMouSh import MoulinShape, TimeStamps, Qin_sinusoidal, Qin_real
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import pandas as pd
import os

secinday = 24*3600
ZERO_KELVIN = 273.15
timestep = 300
supraglacial_baseflow = 0.1

def find_nearest_idx(array, value):
    """Finds the nearest value in an array and outputs a index.
    This function was found in 
    https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array

    Parameters:
    -----------
    array: array to be looked into
    value: single value to look for into the array

    Output:
    -------
    index of the closest value in the array
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def Fill_dict(Q_csv_name,head_csv_name,timestep):
    dictionnary = defaultdict(list)
    
    tmp1 = pd.read_csv(Q_csv_name)
    tmp1 = tmp1.dropna()
    Qin = tmp1.Qm3s.to_numpy() + supraglacial_baseflow
    Qtime = tmp1.SOY.to_numpy()
        
    #time array in seconds
    dictionnary['meltwater_time'] = TimeStamps(Qtime[0],Qtime[-1],timestep)
    dictionnary['meltwater_data'] = Qin_real(dictionnary['meltwater_time'], Qin, Qtime)
    
    tmp2 = pd.read_csv(head_csv_name)
    tmp2 = tmp2.dropna()
    dictionnary['h_real'] = tmp2.head_bed.to_numpy()
    dictionnary['t_real'] = tmp2.soy.to_numpy()
    return dictionnary

def Fill_dict_real(head_csv_name):
    dictionnary = defaultdict(list)
    tmp2 = pd.read_csv(head_csv_name)
    tmp2 = tmp2.dropna()
    dictionnary['h_real'] = tmp2.head_bed.to_numpy()
    dictionnary['t_real'] = tmp2.soy.to_numpy()
    return dictionnary


#temperature profile -- same for all
tmp = pd.read_csv('Field_Data/temperature_foxx1.csv')
temperature_profile = tmp.temperature.to_numpy() #np.array([ZERO_KELVIN, ZERO_KELVIN])#

#field data dictionnaries
tmp = pd.read_csv('Field_Data/surface_discharge_andrews2014_1h_smooth.csv')
tmp = tmp.dropna()

m3 = Fill_dict_real('Field_Data/head_m3.csv')
Qin = tmp.m3_m3s_1h_24hS.to_numpy() + supraglacial_baseflow
Qin[Qin<0]=0
Qtime = tmp.UTC_SOY_1h.to_numpy()
m3['meltwater_time'] = TimeStamps(Qtime[0],Qtime[-1],300)
m3['meltwater_data'] = Qin_real(m3['meltwater_time'], Qin, Qtime)

m4 = Fill_dict_real('Field_Data/head_m4.csv')
Qin = tmp.m4_m3s_1h_24hS.to_numpy() + supraglacial_baseflow
Qin[Qin<0]=0
Qtime = tmp.UTC_SOY_1h.to_numpy()
m4['meltwater_time'] = TimeStamps(Qtime[0],Qtime[-1],300)
m4['meltwater_data'] = Qin_real(m3['meltwater_time'], Qin, Qtime)

foxx = Fill_dict_real('Field_Data/head_mf.csv')
Qin = tmp.mF_m3s_1h_24hS.to_numpy() + supraglacial_baseflow
Qin[Qin<0]=0
Qtime = tmp.UTC_SOY_1h.to_numpy()
foxx['meltwater_time'] = TimeStamps(Qtime[0],Qtime[-1],300)
foxx['meltwater_data'] = Qin_real(m3['meltwater_time'], Qin, Qtime)

jeme = Fill_dict('Field_Data/surface_melt_jeme.csv','Field_Data/head_jeme.csv',timestep)
pira = Fill_dict('Field_Data/surface_melt_pira.csv','Field_Data/head_pira.csv',timestep)
radi = Fill_dict('Field_Data/surface_melt_radi.csv','Field_Data/head_radi.csv',timestep)

#glacier properties (from readme)
m3['ice_thickness']=560
m4['ice_thickness']=540
foxx['ice_thickness']=620
jeme['ice_thickness']=500
pira['ice_thickness']=500
radi['ice_thickness']=700

jeme['baseflow_one'] = np.ones(len(jeme['meltwater_time']))
pira['baseflow_one'] = np.ones(len(pira['meltwater_time']))
radi['baseflow_one'] = np.ones(len(radi['meltwater_time']))
m3['baseflow_one'] = np.ones(len(m3['meltwater_time']))
m4['baseflow_one'] = np.ones(len(m4['meltwater_time']))
foxx['baseflow_one'] = np.ones(len(foxx['meltwater_time']))

# estimate initial subglacial channel area by considering water velocity at 1m/s
# m3['initial_subglacial_area'] = np.mean(m3['meltwater_data'])
# m4['initial_subglacial_area'] = np.mean(m4['meltwater_data'])
# foxx['initial_subglacial_area'] = np.mean(foxx['meltwater_data'])
# jeme['initial_subglacial_area'] = np.mean(jeme['meltwater_data'])
# pira['initial_subglacial_area'] = np.mean(pira['meltwater_data'])
# radi['initial_subglacial_area'] = np.mean(radi['meltwater_data'])

m3['name'] = 'M3'
m4['name'] = 'M4'
foxx['name'] = 'FOXX'
jeme['name'] = 'JEME'
pira['name'] = 'PIRA'
radi['name'] = 'RADI'


m3['time_lim'] = [150,250]
m4['time_lim'] = [150,250]
foxx['time_lim'] = [150,250]
jeme['time_lim'] = [182,250]
pira['time_lim'] = [200,250]
radi['time_lim'] = [208,250]

#surface slope from hoffman 2016 supplemental
regional_surface_slope = 0.01 
channel_length = 25000

del Qin, Qtime, tmp

# # BRUT
# for dataset in [m3,foxx,jeme,radi]:
#     #initiate moulin
#     moulin_sim = MoulinShape(channel_length = channel_length,
#                             temperature_profile = temperature_profile,                   
#                             ice_thickness = dataset['ice_thickness'],
#                             regional_surface_slope = regional_surface_slope,
#                             initial_subglacial_area = dataset['initial_subglacial_area'])
#     portion = (dataset['meltwater_time'] > dataset['time_lim'][0]*secinday) & (dataset['meltwater_time'] < dataset['time_lim'][1]*secinday)
#     time = np.arange(dataset['meltwater_time'][portion][0],dataset['meltwater_time'][portion][-1],timestep*12)
#     #create directory to save figures
#     directory = 'figure_movie_AGU/'+dataset['name']+'_brut'
#     os.mkdir(directory)
    
#     #calculate simulation
#     for idx,t in enumerate(time):
#         meltwater = dataset['meltwater_data'][dataset['meltwater_time']==t]
#         moulin_sim.run1step(t,timestep,meltwater)
        
#     #plot simulation
#     for idx,t_end in enumerate(time):
#     #t_end = time[-1]
#         t_start = t_end-time[-1]+time[0]
#         fig = plt.figure(figsize=(13,5),dpi=150)
#         fig.suptitle(dataset['name'], fontsize=16)
#         moulin_sim.plot_AGU_4(fig,t_start, t_end,
#                              dataset['t_real'],dataset['h_real'],
#                              spine_head_min=200,
#                              ground_depth=-60,
#                              Q_lim = [min(dataset['meltwater_data']),max(dataset['meltwater_data'])],
#                              SC_lim = [min(moulin_sim.dict['subglacial_radius']),max(moulin_sim.dict['subglacial_radius'])],
#                              display_baseflow = False)

#         plt.savefig(directory + '/' + dataset['name'] + '_brut_%d.png'%idx)
#         plt.clf()
#         plt.close(fig)
#     print(dataset['name'],' --> done')
        
#

#subglacial baseflow

#%%
#add meltwater delay
# jeme['meltwater_time'] = m3['meltwater_time']+5*3600

initial_subglacial_area = (np.pi*1**2)/2
fluidity_coefficient_SUB = 6e-24
channel_length = 70000 
creep_factor = 3 
baseflow = 20 
friction = 0.1 
dataset = m3 
path = 'figure_movie_AGU/'
params =  dataset['name']+'_baseflow%d'%baseflow + '_channel%d'%channel_length + '_creep%d'%creep_factor + '_friction%e.1'%friction + 'fluidity_coefficient_SUB%e.1'%fluidity_coefficient_SUB 

#directory = dataset['name'] + params + '/'
#os.mkdir(path + directory)
                                       
moulin_sim = MoulinShape(channel_length = channel_length,                        temperature_profile = temperature_profile,                   
                        ice_thickness = dataset['ice_thickness'],
                        regional_surface_slope = regional_surface_slope,
                        initial_subglacial_area = initial_subglacial_area, 
                        friction_factor_SUB = friction,
                        creep_enhancement_factor = creep_factor,
                        fluidity_coefficient_SUB = fluidity_coefficient_SUB,                              
                        )

portion = (dataset['meltwater_time'] > dataset['time_lim'][0]*secinday) & (dataset['meltwater_time'] < dataset['time_lim'][1]*secinday)
time = dataset['meltwater_time'][portion]
#create directory to save figures



#calculate simulation
for idx,t in enumerate(time):
    meltwater = dataset['meltwater_data'][dataset['meltwater_time']==t]
    moulin_sim.run1step(t,timestep,meltwater,
                        subglacial_baseflow = baseflow, #dataset['meltwater_data'][dataset['meltwater_time']==t]
                        potential_drop=False,
                        open_channel_melt=True,
                        min_radius = 0.1 
                        )
    
#plot simulation
#for idx,t_end in enumerate(time):
t_end = time[-1]
t_start = t_end-time[-1]+time[0]
fig = plt.figure(figsize=(13,5),dpi=150)
fig.suptitle(dataset['name'], fontsize=16)
moulin_sim.plot_AGU_4(fig,t_start, t_end,
                     dataset['t_real'],dataset['h_real'],
                     spine_head_min=200,
                     ground_depth=-60,
                     Q_lim = [min(dataset['meltwater_data']),max(dataset['meltwater_data'])],
                     SC_lim = [min(moulin_sim.dict['subglacial_radius']),max(moulin_sim.dict['subglacial_radius'])],
                     display_baseflow = False)

#plt.savefig(path + directory + '_no%d.png'%idx)
plt.savefig(path + params + '_no%d.png'%idx)
# plt.clf()
# plt.close(fig)
del fig
del moulin_sim



# #%%
# param_csv = pd.read_csv('param_AGU.csv')


# for idx_param in np.arange(len(param_csv)):

#     initial_subglacial_area = (np.pi*param_csv['radius'][idx_param]**2)/2
#     fluidity_coefficient_SUB = param_csv['A'][idx_param] #6e-24
#     channel_length = param_csv['channel_length'][idx_param] #25000 #in [10000,20000,30000,40000,50000,60000]:
#     creep_factor = param_csv['creep_factor'][idx_param] #3 #in [0,3,5,10]:
#     baseflow = param_csv['baseflow'][idx_param] #3 #in [0,1,2,3,4,5,6,7]:
#     friction = param_csv['friction_sub'][idx_param] #0.1 #in [0.001,0.01,0.1,1]:
#     dataset = jeme #in [m3,jeme]:
#     path = 'figure_movie_AGU/'
#     params =  dataset['name']+'_baseflow%d'%baseflow + '_channel%d'%channel_length + '_creep%d'%creep_factor + '_friction%e.1'%friction + 'fluidity_coefficient_SUB%e.1'%fluidity_coefficient_SUB 

#     #directory = dataset['name'] + params + '/'
#     #os.mkdir(path + directory)
                                           
#     moulin_sim = MoulinShape(channel_length = channel_length,
#                             temperature_profile = temperature_profile,                   
#                             ice_thickness = dataset['ice_thickness'],
#                             regional_surface_slope = regional_surface_slope,
#                             initial_subglacial_area = initial_subglacial_area, 
#                             friction_factor_SUB = friction,
#                             creep_enhancement_factor = creep_factor,
#                             fluidity_coefficient_SUB = fluidity_coefficient_SUB                                
#                             )
    
#     portion = (dataset['meltwater_time'] > dataset['time_lim'][0]*secinday) & (dataset['meltwater_time'] < dataset['time_lim'][1]*secinday)
#     time = dataset['meltwater_time'][portion]
#     #create directory to save figures
    
    
    
#     #calculate simulation
#     for idx,t in enumerate(time):
#         meltwater = dataset['meltwater_data'][dataset['meltwater_time']==t]
#         moulin_sim.run1step(t,timestep,meltwater,
#                             subglacial_baseflow = baseflow#dataset['meltwater_data'][dataset['meltwater_time']==t]
#                             )
        
#     #plot simulation
#     #for idx,t_end in enumerate(time):
#     t_end = time[-1]
#     t_start = t_end-time[-1]+time[0]
#     fig = plt.figure(figsize=(13,5),dpi=150)
#     fig.suptitle(dataset['name'], fontsize=16)
#     moulin_sim.plot_AGU_4(fig,t_start, t_end,
#                          dataset['t_real'],dataset['h_real'],
#                          spine_head_min=200,
#                          ground_depth=-60,
#                          Q_lim = [min(dataset['meltwater_data']),max(dataset['meltwater_data'])],
#                          SC_lim = [min(moulin_sim.dict['subglacial_radius']),max(moulin_sim.dict['subglacial_radius'])],
#                          display_baseflow = False)
    
#     #plt.savefig(path + directory + '_no%d.png'%idx)
#     plt.savefig(path + params + '_no%d.png'%idx)
#     # plt.clf()
#     # plt.close(fig)
#     del fig
#     del moulin_sim

