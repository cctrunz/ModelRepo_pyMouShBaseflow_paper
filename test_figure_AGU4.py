from PyMouSh import MoulinShape, TimeStamps, Qin_sinusoidal, Qin_real
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import pandas as pd
import os

secinday = 24*3600
ZERO_KELVIN = 273.15
timestep = 300
supraglacial_baseflow = 0.01


def find_nearest(array, value):
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
    return array[idx]

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
regional_surface_slope = 0.01 
channel_length = 25000

jeme['baseflow_one'] = np.ones(len(jeme['time']))
pira['baseflow_one'] = np.ones(len(pira['time']))
radi['baseflow_one'] = np.ones(len(radi['time']))
m3['baseflow_one'] = np.ones(len(m3['time']))
m4['baseflow_one'] = np.ones(len(m4['time']))
foxx['baseflow_one'] = np.ones(len(foxx['time']))

# estimate initial subglacial channel area by considering water velocity at 1m/s
m3['initial_subglacial_area'] = np.mean(m3['meltwater_input'])
m4['initial_subglacial_area'] = np.mean(m4['meltwater_input'])
foxx['initial_subglacial_area'] = np.mean(foxx['meltwater_input'])
jeme['initial_subglacial_area'] = np.mean(jeme['meltwater_input'])
pira['initial_subglacial_area'] = np.mean(pira['meltwater_input'])
radi['initial_subglacial_area'] = np.mean(radi['meltwater_input'])

m3['name'] = 'Moulin 3'
m4['name'] = 'Moulin 4'
foxx['name'] = 'Foxx moulin'
jeme['name'] = 'Jeme moulin'
pira['name'] = 'Pirate moulin'
radi['name'] = 'Radical moulin'

del Qin, Qtime, tmp
##idealized Qin for test purposes
# Qin_mean = 1
# dQ = 0.1 
# meltwater_input = Qin_sinusoidal(time,Qin_mean, dQ)


for dataset in [m3,foxx,jeme,radi]:
    moulin_sim = MoulinShape(channel_length = channel_length,
                            temperature_profile = temperature_profile,                   
                            ice_thickness = dataset['ice_thickness'],
                            regional_surface_slope = regional_surface_slope,
                            initial_subglacial_area = dataset['initial_subglacial_area'])
    start = 0
    end = 3
    time = dataset['time'][int(start*secinday/300):int(end*secinday/300)]
    
    for idx,t in enumerate(time):
        meltwater = dataset['meltwater_input'][idx]
        moulin_sim.run1step(t,timestep,meltwater)
    
    
    
    #for idx,t_end in enumerate(time):
    t_end = time[-1]
    t_start = t_end-time[-1]+time[0]
    fig = plt.figure(figsize=(13,5),dpi=150)
    fig.suptitle(dataset['name'], fontsize=16)
    moulin_sim.plot_AGU_4(fig,t_start, t_end,
                         dataset['t_real'],dataset['h_real'],
                         spine_head_min=200,
                         ground_depth=-60,
                         Q_lim = [min(dataset['meltwater_input']),max(dataset['meltwater_input'])],
                         SC_lim = [min(moulin_sim.dict['subglacial_radius']),max(moulin_sim.dict['subglacial_radius'])],
                         display_baseflow = False)
    os.mkdir('figure_movie_AGU/'+dataset['name']+'_brut')
    plt.savefig('jeme_brut_%d.png'%idx)

