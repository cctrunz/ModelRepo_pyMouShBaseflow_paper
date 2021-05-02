from pyMouSh import MoulinShape, TimeStamps, Qin_constant, Qin_sinusoidal, Qin_real, calculate_h_S_schoof, find_nearest
import numpy as np
from datetime import datetime as dt  
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
from matplotlib.patches import ConnectionPatch
import pandas as pd
import pickle
import seaborn as sns




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
#stage_data = stage.stage.to_numpy()
#stage_doy = stage.DOY.to_numpy()

#import JEME dye
dye = pd.read_csv('Field_Data/JEME_Qdye.csv', index_col=0, parse_dates=True)
elapsed = dye.index-dt(2017, 1, 1,00,00,00)
dye['SOY']=elapsed.total_seconds()
dye['DOY']=dye.SOY/3600/24  
#Qdye = dye['Discharge (m^3/s)'].to_numpy()
#dye_soy = dye.SOY.to_numpy()
#dye_doy = dye.DOY.to_numpy()

#import melt model inputs
weather = pd.read_csv('Field_Data/LOWC17_MELT_NEW.csv', index_col=0, parse_dates=True)
elapsed = weather.index-dt(2017, 1, 1,00,00,00)
weather['SOY']=elapsed.total_seconds()
weather['DOY']=weather.SOY/3600/24  
#wea_soy = weather.SOY.to_numpy()
#wea_doy = weather.DOY.to_numpy()
#air_temp = weather.Temp.to_numpy()
#albedo = weather.albedo.to_numpy()
#sw_down = weather.SW_down.to_numpy()

#Import meltwater input calculated with meltmodel for JEME
#drainage_basin = 0.9e6 #in m^2
meltmodel = pd.read_csv('Field_Data/jeme_hourly_melt_and_qsuh.csv')
meltmodel = meltmodel.dropna()
#melt = meltmodel.melt_mmperh_smooth_6h.to_numpy()# mm/h
#Qmelt = (melt/1e3/3600) * drainage_basin
#qsuh = meltmodel.q_suh_mmh #mm/h
#Qsuh = (qsuh/1e3/3600) * drainage_basin
#meltmodel_soy = meltmodel.SOY.to_numpy() #add routing delay -- Jess found 2h not 4. investigate why
#meltmodel_doy = meltmodel.DOY.to_numpy()

#Import head measurement for JEME
jeme_moulin = pd.read_csv('Field_Data/head_jeme.csv')
jeme_moulin = jeme_moulin.dropna()
h_real = jeme_moulin.head_bed.to_numpy()
h_soy = jeme_moulin.soy.to_numpy()
h_doy = h_soy/secinday

#set timeserie length and timestep for interpolation
time_start = meltmodel.SOY[int(2*secinday/900)]  
time_end = time_start + 65*secinday
timestep = 10*60#300 #seconds to reduce number of saved data
time = TimeStamps(time_start,time_end,timestep)

#calculate meltwater input (interpolate)
meltwater_input = Qin_real(time, meltmodel.Qm3s.to_numpy()+0.01,  meltmodel.SOY.to_numpy())

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


#%%

fig= plt.figure(figsize=(4,5))
# grid = plt.GridSpec(8,1, hspace=-5)
# ax1 = fig.add_subplot(grid[0, 0]) 
# ax2 = fig.add_subplot(grid[1, 0])#,sharex=ax1)  
# ax3 = fig.add_subplot(grid[2, 0])#,sharex=ax1)  
# ax4 = fig.add_subplot(grid[3, 0])#,sharex=ax1)  
# ax5 = fig.add_subplot(grid[4:8, 0])#,sharex=ax1)  

grid = plt.GridSpec(10,1)
ax1 = fig.add_subplot(grid[0:2, 0]) 
ax2 = fig.add_subplot(grid[1:3, 0]) 
ax3 = fig.add_subplot(grid[2:4, 0]) 
ax4 = fig.add_subplot(grid[3:6, 0]) 
ax5 = fig.add_subplot(grid[5:9, 0])
     
mean = 0.15
App = 0.22
xlim = [200,210]

#idealized meltwater input arrays
sin_melt_input = Qin_sinusoidal(time, mean, App/2, 15)

#Weather
sns.lineplot(x='DOY', y='Temp', data=weather,ax=ax1, color='black')
ax1.fill_between(weather.DOY.to_numpy(),0,weather.Temp.to_numpy(),color='lightgrey')
ax1.set_ylim(-5,10)
ax1.set_yticks([0,4,8])
ax1.spines['left'].set_bounds([0,8])

sns.lineplot(x='DOY', y='albedo', data=weather,ax=ax2, color='grey')
ax2.set_ylim(0.2,0.8)
ax2.set_yticks([0.7,0.5,0.3])
ax2.spines['right'].set_bounds([0.7,0.3])

sns.lineplot(x='DOY', y='SW_down', data=weather,ax=ax3, color='black')
ax3.set_ylim(-50,1000)
ax3.set_yticks([0,300,600,900])
ax3.spines['left'].set_bounds([0,900])
#stage
#sns.lineplot(x='DOY', y='stage', data=stage,ax=ax4, color='grey')
ax4.plot(stage['2017-07-22':].DOY,stage.stage['2017-07-22':], color='grey')
ax4.set_ylabel('Stream height')
ax4.set_ylim(2.5,3.5)
ax4.set_yticks([2.5,2.7,2.9,3.1])
ax4.spines['right'].set_bounds([2.5,3.1])


#move axis to right

sns.despine(ax=ax1 , bottom=True, right=True)
sns.despine(ax=ax2, left=True, bottom=True, right=False)
sns.despine(ax=ax3, bottom=True)
sns.despine(ax=ax4, left=True, bottom=True, right=False)
sns.despine(ax=ax5, offset=(10,0))

ax1.tick_params(bottom=False, labelbottom=False)
ax2.tick_params(bottom=False, labelbottom=False, color='grey') 
ax3.tick_params(bottom=False, labelbottom=False) 
ax4.tick_params(bottom=False, labelbottom=False, color='grey')    

ax1.set(xlabel=None)
ax2.set(xlabel=None)
ax3.set(xlabel=None)
ax4.set(xlabel=None)
ax5.set(xlabel=None)

ax2.yaxis.set_label_position("right")
ax4.yaxis.set_label_position("right")

ax2.spines['right'].set_color('grey')
ax4.spines['right'].set_color('grey')

#rectangle amplitude of oscillation
ax5.axhspan(mean-App/2,mean+App/2, color='gray', alpha=0.15)
ax5.axhline(y=mean, color='black', linestyle='--')
ax5.add_artist(ConnectionPatch((200,mean-App/2),(200,mean+App/2),'data','data', shrinkB=0,
                arrowstyle="<->", color='black', linewidth=2))
ax5.text(199.8,0.28,'$A_{in}$',fontsize=10)
ax5.text(210.1,0.13,'$\overline{Q}_{in}$',fontsize=10)

#Melt model
sns.lineplot(data=meltmodel, x="DOY", y="Qm3s", ax=ax5, color=c_input)
#Measured
ax5.errorbar( dye.DOY.to_numpy(), dye['Discharge (m^3/s)'].to_numpy(),yerr=dye['Discharge (m^3/s)'].to_numpy()*0.25, 
            label='Measured $Q_{in}$', color = 'orangered', linewidth=1, ecolor='darksalmon',zorder=1)
#ax.errorbar(meltmodel_doy, Qsuh, yerr=Qsuh*0.25, 
#            color=c_input, label='Modeled $Q_{in}$',linewidth=2, ecolor=c_input_light, zorder=0)#'steelblue','lightblue'

#quartiles3 = dye['Discharge (m^3/s)'].to_numpy() + dye['Discharge (m^3/s)'].to_numpy() *0.25
#quartiles1 = dye['Discharge (m^3/s)'].to_numpy() - dye['Discharge (m^3/s)'].to_numpy() *0.25
#sns.lineplot(data=dye, x="DOY", y="Discharge (m^3/s)", ax=ax5)
#ax5.fill_between(dye['Discharge (m^3/s)'].to_numpy(), quartiles1, quartiles3, alpha=0.3, color='orangered'); 

ax1.set_xlim(xlim)
ax2.set_xlim(xlim)
ax3.set_xlim(xlim)
ax4.set_xlim(xlim)
ax5.set_xlim(xlim)


ax5.set_ylabel('$Q_{in}$ ($m^3s^{-1}$)')
ax5.set_xlabel('Days of year 2017')
ax5.set_ylim(0,0.5)

ax1.patch.set_alpha(0)
ax2.patch.set_alpha(0)
ax3.patch.set_alpha(0)
ax4.patch.set_alpha(0)
ax5.patch.set_alpha(0)


#legend
elements = [Line2D([0], [0], color='orangered', ls='-', lw=2, label='Measured'),
            Line2D([0], [0], color=c_input, ls='-', lw=2, label='Modeled'),
            Line2D([0], [0], color='black', lw=1, ls=':', label='Idealized'),]
ax5.legend(handles=elements, loc=2, bbox_to_anchor=(-0.05, 1), labelspacing=0,  prop={'size': 8})



plt.savefig('Figure_TC/stream_discharge.pdf',bbox_inches = 'tight')







