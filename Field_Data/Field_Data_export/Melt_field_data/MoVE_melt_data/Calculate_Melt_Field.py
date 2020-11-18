

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


from utils_jessica_oct2020 import c_rolling
from datetime import datetime as dt


jeme_basin_area = 8e5 #m2
pira_basin_area = 8e5 #m2
radi_basin_area = 6e6 #m2


def create_csv(df,basin_area,year,filename='filename.csv'):
    '''smooth data with a rolling window'''
    # center a 6h rolling mean
    df['melt_rate_smooth_6h'] = c_rolling(df.melt_rate, '6H')
    '''transform units from mm/h to m3/s'''
    df['Qm3s'] = df.melt_rate_smooth_6h * basin_area /1e3/3600         
    '''add DOY and SOY'''
    elapsed = df.index-dt(year, 1, 1,00,00,00)
    df['SOY']=elapsed.total_seconds()
    df['DOY']=df.SOY/3600/24    
    '''export to csv'''
    df.to_csv(filename)
    
create_csv(lowc17_melt,jeme_basin_area,2017,filename='surface_melt_jeme.csv')
create_csv(lowc18_melt,pira_basin_area,2017,filename='surface_melt_pira.csv')
create_csv(high17_melt,radi_basin_area,2017,filename='surface_melt_radi.csv')


