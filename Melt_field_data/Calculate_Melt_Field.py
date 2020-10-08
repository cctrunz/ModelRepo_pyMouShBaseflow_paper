

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

#%% plot

import matplotlib.pyplot as plt

plt.figure()
plt.plot(lowc17_melt)
plt.plot(lowc18_melt)
plt.plot(high17_melt)
plt.plot(high18_melt)



#%%
#smooth data with a rolling window.

lowc17_melt_rolling = lowc17_melt.rolling(10).mean()
lowc18_melt_rolling = lowc18_melt.rolling(10).mean()
high17_melt_rolling = high17_melt.rolling(10).mean()
high18_melt_rolling = high18_melt.rolling(10).mean()

plt.figure()
plt.plot(lowc17_melt,color='blue')
plt.plot(lowc18_melt,color='blue')
plt.plot(high17_melt,color='blue')
plt.plot(high18_melt,color='blue')
plt.plot(lowc17_melt_rolling, color='red')
plt.plot(lowc18_melt_rolling, color='red')
plt.plot(high17_melt_rolling, color='red')
plt.plot(high18_melt_rolling, color='red')

#%%
#save to csv file

lowc17_melt_rolling.to_csv('lowc17_melt_rolling.csv')
lowc18_melt_rolling.to_csv('lowc18_melt_rolling.csv')
high17_melt_rolling.to_csv('high17_melt_rolling.csv')
high18_melt_rolling.to_csv('high18_melt_rolling.csv')
