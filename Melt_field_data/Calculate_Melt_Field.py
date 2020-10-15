

'''Initial code and functions created by Jessica Mejia
Modified 2020 by Celia Trunz'''

'''Read in concatenated weather data csv files (HOBO)
create object for each set of weather data
determine incoming vs reflected solar radiation and
name accordingly. 
Calculate melt rate (mm w.e. h^-1) using the model of 
Pellicciotti et al., 2015 using daily albedo values
if avaiable.
Where incoming and reflected solar radiation sensor data
is not available use incoming data from other station if
available. If no reflected radiation collected locally, use
a standard albedo of 0.7 for entire timeseries'''


from melt_model_jessica_Oct2020 import weather_data

lowc17= weather_data('LOWC17_WEA_CAT.csv')
lowc17_melt  = weather_data.calc_melt(lowc17)

high17 = weather_data('HIGH17_WEA_CAT.csv')
high17_melt  = weather_data.calc_melt(high17)


# 2018 weather data
lowc18 = weather_data('LOWC18_WEA_CAT.csv')
lowc18_melt = weather_data.calc_melt(lowc18)

# 2018 high camp data is missing, use solar rad from low camp
high18 = weather_data('HIGH18_WEA_CAT.csv')
high18_melt = high18.calc_melt(incoming_shortwave_radiation=lowc18.solar_corrected)


#%% 
'''plot melt rate'''

import matplotlib.pyplot as plt

plt.figure()
plt.plot(lowc17_melt,label='LC17')
plt.plot(lowc18_melt,label='LC18')
plt.plot(high17_melt,label='HC17')
plt.plot(high18_melt,label='HC18')
plt.ylabel('melt w.e. ($mm/h$)')
plt.legend()



#%% 
'''Transform to discharge
drainage areas
--------------
jeme 2017 : 1e6 m2 or 3e5
pira 2018 : 1e6 m2 or 3e5

radi 2017 : 6e6 m2

'''

jeme_basin_area = 8e5 #m2
pira_basin_area = 8e5 #m2
radi_basin_area = 6e6 #m2


#transform units from mm/h to m3/s

lowc17_Q_jeme = lowc17_melt * jeme_basin_area /1e3/3600 
lowc18_Q_pira = lowc18_melt * pira_basin_area /1e3/3600 
high17_Q_radi = high17_melt * radi_basin_area /1e3/3600 
high18_Q_radi = high18_melt * radi_basin_area /1e3/3600 

'''Import pira discharge data to compare'''
from pandas import read_csv

pira_ir=read_csv('PIRA18_IR.csv',index_col=0,parse_dates=True)
pira_dye=read_csv('PIRA_Qdye.csv',index_col=0,parse_dates=True)



plt.figure()
plt.plot(lowc17_Q_jeme,label='LC17')
plt.plot(lowc18_Q_pira,label='LC18')
plt.plot(high17_Q_radi,label='HC17')
plt.plot(high18_Q_radi,label='HC18')
pira_ir['Discharge smoothed (m^3/s)'].plot(label='IceRay')
pira_dye['Discharge (m^3/s)'].plot(style='.',label='Dye')
plt.ylabel('Discharge ($m^3/s$)')
plt.legend()



#%%


#%%
'''smooth data with a rolling window'''


lowc17_Q_jeme_rolling = lowc17_Q_jeme.rolling(15,min_periods=1).mean() 
lowc18_Q_pira_rolling = lowc18_Q_pira.rolling(15,min_periods=1).mean() 
high17_Q_radi_rolling = high17_Q_radi.rolling(15,min_periods=1).mean()
high18_Q_radi_rolling = high18_Q_radi.rolling(15,min_periods=1).mean() 

plt.figure()
plt.plot(lowc17_Q_jeme,color='blue')
plt.plot(lowc18_Q_pira,color='blue')
plt.plot(high17_Q_radi,color='blue')
plt.plot(high18_Q_radi,color='blue')
plt.plot(lowc17_Q_jeme_rolling, color='red')
plt.plot(lowc18_Q_pira_rolling, color='red')
plt.plot(high17_Q_radi_rolling, color='red')
plt.plot(high18_Q_radi_rolling, color='red')


'''create elapsed time array'''
elapsed_l17 = lowc17_Q_jeme_rolling.index-lowc17_Q_jeme_rolling.index[0]
elapsed_l18 = lowc18_Q_pira_rolling.index-lowc18_Q_pira_rolling.index[0]
elapsed_h17 = high17_Q_radi_rolling.index-high17_Q_radi_rolling.index[0]
elapsed_h18 = high18_Q_radi_rolling.index-high18_Q_radi_rolling.index[0]
lowc17_Q_jeme_rolling['Seconds']=elapsed_l17.total_seconds()
lowc18_Q_pira_rolling['Seconds']=elapsed_l18.total_seconds()
high17_Q_radi_rolling['Seconds']=elapsed_h17.total_seconds()
high18_Q_radi_rolling['Seconds']=elapsed_h18.total_seconds()


#%%
#save to csv file

lowc17_Q_jeme_rolling.to_csv('lowc17_Q_jeme_rolling.csv')
lowc18_Q_pira_rolling.to_csv('lowc18_Q_pira_rolling.csv')
high17_Q_radi_rolling.to_csv('high17_Q_radi_rolling.csv')
high18_Q_radi_rolling.to_csv('high18_Q_radi_rolling.csv')
