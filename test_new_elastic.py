#%pylab inline
from PyMouSh import MoulinShape, TimeStamps, Qin_constant, Qin_sinusoidal, Qin_real, calculate_h_S_schoof, find_nearest
import numpy as np
from datetime import datetime as dt  
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
from matplotlib.patches import ConnectionPatch
import pandas as pd
import pickle


### 
secinday = 24*3600
ZERO_KELVIN = 273.15
#paramters

#temperature profile based of foxx data -- Lauren
tmp = pd.read_csv('Field_Data/temperature_foxx1.csv')
temperature_profile = tmp.temperature.to_numpy()#np.array([ZERO_KELVIN, ZERO_KELVIN])

moulin_radii = 0.2
initial_subglacial_area = 1 
regional_surface_slope = 0.01
channel_length = 25e3
ice_thickness = 500
initial_head = ice_thickness

#import JEME stage
stage = pd.read_csv('Field_Data/LOWC17_STR_1.csv', index_col=0, parse_dates=True)
elapsed = stage.index-dt(2017, 1, 1,00,00,00)
stage['SOY']=elapsed.total_seconds()
stage['DOY']=stage.SOY/3600/24  

#import JEME dye
dye = pd.read_csv('Field_Data/JEME_Qdye.csv', index_col=0, parse_dates=True)


#import melt model inputs
weather = pd.read_csv('Field_Data/LOWC17_MELT_NEW.csv', index_col=0, parse_dates=True)
elapsed = weather.index-dt(2017, 1, 1,00,00,00)
weather['SOY']=elapsed.total_seconds()
weather['DOY']=weather.SOY/3600/24  

#Import meltwater input calculated with meltmodel for JEME
meltmodel = pd.read_csv('Field_Data/JEME_QSUH.csv', index_col=0, parse_dates=True)
meltmodel = meltmodel.dropna()
meltmodel = meltmodel['2017/05/30':'2017/09/01']#['2017/07/19':'2017/09/01']#['2017/05/30':'2017/09/01']


#Import head measurement for JEME
jeme_moulin = pd.read_csv('Field_Data/head_jeme.csv')
jeme_moulin = jeme_moulin.dropna()
h_real = jeme_moulin.head_bed.to_numpy()
h_soy = jeme_moulin.soy.to_numpy()
h_doy = h_soy/secinday

#set timeserie length and timestep for interpolation
time_start = meltmodel.soy[0]#meltmodel.soy[int(2*secinday/900)]  
time_end = meltmodel.soy[-1]#time_start + 65*secinday
timestep = 10*60#300 #seconds to reduce number of saved data
time = TimeStamps(time_start,time_end,timestep)

#calculate meltwater input (interpolate)
meltwater_input = Qin_real(time, meltmodel.Q.to_numpy()+0.02,  meltmodel.soy.to_numpy())

#calculate baseflow
bf_mean = 0.5
bf_amp = 0.05 
shift0h = 0.42 * secinday
shift12h = 0 * secinday # sin is already shifted from input
baseflow_shift0 = Qin_sinusoidal(time,bf_mean, bf_amp, shift=shift0h)
baseflow_shift12 = Qin_sinusoidal(time,bf_mean, bf_amp, shift=shift12h)

c_baserock = 'peru'#'tan'#'sienna'##'grey'
c_cyl10 = (0.00392156862745098, 0.45098039215686275, 0.6980392156862745)
c_cyl1 = 'lightblue'
c_input = 'darkorchid' #'#af8dc3'#(0.8, 0.47058823529411764, 0.7372549019607844)
c_input_light = 'plum'#'#af8dc3'
c_fix = '#01665e'
c_osc = '#35978f'
c_lag = '#80cdc1'
c_2 = '#bf812d'
c_0 = 'grey'
ls_0 = (0, (3, 1, 1, 1, 1, 1))
ls_hreal = '--'
ls_input = (0, (1, 1))

bf0_fix = MoulinShape(moulin_radii = moulin_radii,
                      ice_thickness = ice_thickness,
                      initial_head = initial_head,
                      initial_subglacial_area = np.pi*0.2**2,
                      channel_length = channel_length,
                      regional_surface_slope = regional_surface_slope,
                      temperature_profile = temperature_profile)

for idx,t in enumerate(time) :
    bf0_fix.run1step(t,
                    timestep,
                    meltwater_input[idx],
                     overflow=True,
                    subglacial_baseflow = 0)
    
#%%

plot_data(bf0_fix,h_soy,h_real)
plt.savefig('bf0_fix.pdf')