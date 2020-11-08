# calc melt for foxx weather station (2012)
import pandas as pd
import matplotlib.pyplot as plt
from datatools.melt_model import WeatherStation
# import laurens foxx data
foxx = pd.read_csv("weather_foxx_2012.csv",na_values=[-7999])

foxx['Date'] = foxx['UTC_DOY'].apply(
    lambda x: pd.Timestamp('2012-01-01')+pd.Timedelta(days=x-1))
foxx.set_index('Date',inplace=True)
foxx.index = foxx.index.round(freq='5T')
foxx.rename(columns={'refl_sw_wattsperm2':'Reflected',
                    'incom_sw_wattsperm2':'Solar',
                    'airtemp_C':'Temp'},inplace=True)    
FOXX = WeatherStation(foxx,name='FOXX')
                                     
foxx_melt = FOXX.calc_melt()
plt.figure()
plt.plot(foxx_melt)
plt.ylabel('mm m.w. equiv per hour')     
plt.title('FOXX weather station melt rate')
#%%
foxx_discharge = (foxx_melt/3600/1000)*3028000
plt.figure(figsize=(10,2))
plt.plot(foxx_discharge)
plt.ylabel('$m^3/s$')     
plt.title('FOXX weather station melt rate')
plt.savefig('melt_m3_jessica.png')