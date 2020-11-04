from PyMouSh import MoulinShape, TimeStamps, Qin_sinusoidal, Qin_real
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import pandas as pd

secinday = 24*3600
ZERO_KELVIN = 273.15
timestep = 300
supraglacial_baseflow = 0.01

def Fill_dict(Q_csv_name,head_csv_name,timestep):
    dictionnary = defaultdict(list)
    
    tmp1 = pd.read_csv(Q_csv_name)
    tmp1 = tmp1.dropna()
    Qin = tmp1.Qm3s.to_numpy() + supraglacial_baseflow
    Qtime = tmp1.SOY.to_numpy()
        
    #time array in seconds
    dictionnary['time'] = TimeStamps(Qtime[0],Qtime[-1],timestep)
    dictionnary['meltwater_input'] = Qin_real(dictionnary['time'], Qin, Qtime)
    
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
m3['time'] = TimeStamps(Qtime[0],Qtime[-1],300)
m3['meltwater_input'] = Qin_real(m3['time'], Qin, Qtime)

m4 = Fill_dict_real('Field_Data/head_m4.csv')
Qin = tmp.m4_m3s_1h_24hS.to_numpy() + supraglacial_baseflow
Qin[Qin<0]=0
Qtime = tmp.UTC_SOY_1h.to_numpy()
m4['time'] = TimeStamps(Qtime[0],Qtime[-1],300)
m4['meltwater_input'] = Qin_real(m3['time'], Qin, Qtime)

foxx = Fill_dict_real('Field_Data/head_mf.csv')
Qin = tmp.mF_m3s_1h_24hS.to_numpy() + supraglacial_baseflow
Qin[Qin<0]=0
Qtime = tmp.UTC_SOY_1h.to_numpy()
foxx['time'] = TimeStamps(Qtime[0],Qtime[-1],300)
foxx['meltwater_input'] = Qin_real(m3['time'], Qin, Qtime)

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

#surface slope from hoffman 2016 supplemental
regional_surface_slope = 0.02 
channel_length = 25000

jeme['baseflow_one'] = np.ones(len(jeme['time']))
pira['baseflow_one'] = np.ones(len(pira['time']))
radi['baseflow_one'] = np.ones(len(radi['time']))
m3['baseflow_one'] = np.ones(len(m3['time']))
m4['baseflow_one'] = np.ones(len(m4['time']))
foxx['baseflow_one'] = np.ones(len(foxx['time']))

del Qin, Qtime, tmp
##idealized Qin for test purposes
# Qin_mean = 1
# dQ = 0.1 
# meltwater_input = Qin_sinusoidal(time,Qin_mean, dQ)



#jeme
#%%#######################################################################################
jeme_brut = MoulinShape(channel_length = channel_length,
                        temperature_profile = temperature_profile,                   
                        ice_thickness = jeme['ice_thickness'],
                        regional_surface_slope = regional_surface_slope,
                        initial_subglacial_area = (np.pi*0.2**2)/2)
start = 0
end = 100
time = jeme['time'][int(start*secinday/300):int(end*secinday/300)]

for idx,t in enumerate(time):
    meltwater = jeme['meltwater_input'][idx]
    jeme_brut.run1step(t,timestep,meltwater)



for idx,t_end in enumerate(time):
    t_start = t_end-time[-1]+time[0]
    fig = plt.figure(figsize=(9,6),dpi=150)
    fig.suptitle('JEME (Low camp 2017)', fontsize=16)
    jeme_brut.plot_AGU_4(fig,t_start, t_end,
                         jeme['t_real'],jeme['h_real'],
                         spine_head_min=500,
                         ground_depth=-100,
                         Q_lim = [0,4],
                         SC_lim = [0,0.5],
                         display_baseflow = True)
    plt.savefig('figure_movie_AGU/jeme_brut_%d.png'%idx)



#radi
########################################################################################
radi_brut = MoulinShape(channel_length = channel_length,
                        temperature_profile = temperature_profile,                   
                        ice_thickness = foxx['ice_thickness'],
                        regional_surface_slope = regional_surface_slope,
                        initial_subglacial_area = (np.pi*0.7**2)/2)

start = 0
end = 100
time = radi['time'][int(start*secinday/300):int(end*secinday/300)]

for idx,t in enumerate(time):
    meltwater = radi['meltwater_input'][idx]
    radi_brut.run1step(t,timestep,meltwater)
    
for idx,t_end in enumerate(time):
    t_start = t_end-time[-1]+time[0]
    fig = plt.figure(figsize=(9,6),dpi=150)
    fig.suptitle('RADI (High camp 2017)', fontsize=16)
    radi_brut.plot_AGU_4(fig,t_start, t_end,
                         radi['t_real'],radi['h_real'],
                         spine_head_min=500,
                         ground_depth=-100,
                         Q_lim = [0,4],
                         SC_lim = [0,0.5],
                         display_baseflow = True)
plt.savefig('figure_movie_AGU/radi_brut_%d.png'%idx)

#moulin 3
########################################################################################
m3_brut = MoulinShape(channel_length = channel_length,
                        temperature_profile = temperature_profile,                   
                        ice_thickness = foxx['ice_thickness'],
                        regional_surface_slope = regional_surface_slope,
                        initial_subglacial_area = (np.pi*0.7**2)/2)
start = 180
end = 250
time = m3['time'][int(start*secinday/300):int(end*secinday/300)]

for idx,t in enumerate(time):
    meltwater = m3['meltwater_input'][idx]
    m3_brut.run1step(t,timestep,meltwater)
       
for idx,t_end in enumerate(time):
    t_start = t_end-time[-1]+time[0]
    fig = plt.figure(figsize=(9,6),dpi=150)
    fig.suptitle('Moulin 3', fontsize=16)
    m3_brut.plot_AGU_4(fig,t_start, t_end,
                         m3['t_real'],m3['h_real'],
                         spine_head_min=500,
                         ground_depth=-100,
                         Q_lim = [0,4],
                         SC_lim = [0,0.5],
                         display_baseflow = True)
plt.savefig('figure_movie_AGU/m3_brut_%d.png'%idx)



#moulin foxx
########################################################################################
foxx_brut = MoulinShape(channel_length = channel_length,
                        temperature_profile = temperature_profile,                   
                        ice_thickness = foxx['ice_thickness'],
                        regional_surface_slope = regional_surface_slope,
                        initial_subglacial_area = (np.pi*0.7**2)/2)
start = 180
end = 250
time = foxx['time'][int(start*secinday/300):int(end*secinday/300)]

for idx,t in enumerate(time):
    meltwater = foxx['meltwater_input'][idx]
    foxx_brut.run1step(t,timestep,meltwater)
       
for idx,t_end in enumerate(time):
    t_start = t_end-time[-1]+time[0]
    fig = plt.figure(figsize=(9,6),dpi=150)
    fig.suptitle('Moulin Foxx', fontsize=16)
    foxx_brut.plot_AGU_4(fig,t_start, t_end,
                         foxx['t_real'],foxx['h_real'],
                         spine_head_min=500,
                         ground_depth=-100,
                         Q_lim = [0,4],
                         SC_lim = [0,0.5],
                         display_baseflow = True)
plt.savefig('figure_movie_AGU/foxx_brut_%d.png'%idx)






########################################################################################
########################################################################################
#%% BASEFLOW
########################################################################################
########################################################################################
########################################################################################

#jeme
########################################################################################
jeme_baseflow = MoulinShape(channel_length = channel_length,
                        temperature_profile = temperature_profile,                   
                        ice_thickness = jeme['ice_thickness'],
                        regional_surface_slope = regional_surface_slope,
                        initial_subglacial_area = (np.pi*0.2**2)/2)
start = 0
end = 100
time = jeme['time'][int(start*secinday/300):int(end*secinday/300)]

for idx,t in enumerate(time):
    meltwater = jeme['meltwater_input'][idx]
    jeme_baseflow.run1step(t,timestep,meltwater,
                       subglacial_baseflow=3)
      
for idx,t_end in enumerate(time):
    t_start = t_end-time[-1]+time[0]
    fig = plt.figure(figsize=(9,6),dpi=150)
    fig.suptitle('JEME (Low camp 2017)', fontsize=16)
    jeme_baseflow.plot_AGU_4(fig,t_start, t_end,
                         jeme['t_real'],jeme['h_real'],
                         spine_head_min=500,
                         ground_depth=-100,
                         Q_lim = [0,4],
                         SC_lim = [0,0.5],
                         display_baseflow = True)
plt.savefig('figure_movie_AGU/jeme_brut_baseflow3_%d.png'%idx)

#radi
########################################################################################
radi_baseflow = MoulinShape(channel_length = channel_length,
                        temperature_profile = temperature_profile,                   
                        ice_thickness = foxx['ice_thickness'],
                        regional_surface_slope = regional_surface_slope,
                        initial_subglacial_area = (np.pi*0.7**2)/2)

start = 0
end = 100
time = radi['time'][int(start*secinday/300):int(end*secinday/300)]

for idx,t in enumerate(time):
    meltwater = radi['meltwater_input'][idx]
    radi_baseflow.run1step(t,timestep,meltwater,
                       subglacial_baseflow=3)
    
for idx,t_end in enumerate(time):
    t_start = t_end-time[-1]+time[0]
    fig = plt.figure(figsize=(9,6),dpi=150)
    fig.suptitle('RADI (High camp 2017)', fontsize=16)
    radi_baseflow.plot_AGU_4(fig,t_start, t_end,
                         radi['t_real'],radi['h_real'],
                         spine_head_min=500,
                         ground_depth=-100,
                         Q_lim = [0,4],
                         SC_lim = [0,0.5],
                         display_baseflow = True)
plt.savefig('figure_movie_AGU/radi_brut_baseflow3_%d.png'%idx)

#moulin 3
########################################################################################
m3_baseflow = MoulinShape(channel_length = channel_length,
                        temperature_profile = temperature_profile,                   
                        ice_thickness = foxx['ice_thickness'],
                        regional_surface_slope = regional_surface_slope,
                        initial_subglacial_area = (np.pi*0.7**2)/2)
start = 180
end = 250
time = m3['time'][int(start*secinday/300):int(end*secinday/300)]

for idx,t in enumerate(time):
    meltwater = m3['meltwater_input'][idx]
    m3_baseflow.run1step(t,timestep,meltwater,
                       subglacial_baseflow=3)
       
for idx,t_end in enumerate(time):
    t_start = t_end-time[-1]+time[0]
    fig = plt.figure(figsize=(9,6),dpi=150)
    fig.suptitle('Moulin 3', fontsize=16)
    m3_baseflow.plot_AGU_4(fig,t_start, t_end,
                         m3['t_real'],m3['h_real'],
                         spine_head_min=500,
                         ground_depth=-100,
                         Q_lim = [0,4],
                         SC_lim = [0,0.5],
                         display_baseflow = True)
plt.savefig('figure_movie_AGU/m3_brut_baseflow3_%d.png'%idx)


#moulin foxx
########################################################################################
foxx_baseflow = MoulinShape(channel_length = channel_length,
                        temperature_profile = temperature_profile,                   
                        ice_thickness = foxx['ice_thickness'],
                        regional_surface_slope = regional_surface_slope,
                        initial_subglacial_area = (np.pi*0.7**2)/2)
start = 180
end = 250
time = foxx['time'][int(start*secinday/300):int(end*secinday/300)]

for idx,t in enumerate(time):
    meltwater = foxx['meltwater_input'][idx]
    foxx_baseflow.run1step(t,timestep,meltwater,
                       subglacial_baseflow=3)
       
for idx,t_end in enumerate(time):
    t_start = t_end-time[-1]+time[0]
    fig = plt.figure(figsize=(9,6),dpi=150)
    fig.suptitle('Moulin Foxx', fontsize=16)
    foxx_baseflow.plot_AGU_4(fig,t_start, t_end,
                         foxx['t_real'],foxx['h_real'],
                         spine_head_min=500,
                         ground_depth=-100,
                         Q_lim = [0,4],
                         SC_lim = [0,0.5],
                         display_baseflow = True)
plt.savefig('figure_movie_AGU/foxx_brut_baseflow3_%d.png'%idx)
