import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.dates as mdates

# smooth1h = pd.read_csv('surface_discharge_andrews2014_1h_smooth.csv')
# six_hours = pd.read_csv('surface_discharge_andrews2014_6h.csv')

melt = pd.read_csv('surface_melt_pira.csv',parse_dates=['Date'])
dye = pd.read_csv('PIRA_Qdye.csv',parse_dates=['Time (UTC-0)'])
iceray = pd.read_csv('PIRA18_IR.csv',parse_dates=['Time (UTC-0)'])

#%%
fig = plt.figure(figsize=(6,2))
ax = fig.add_subplot(111)
plt.plot(dye['Time (UTC-0)'],dye['Discharge (m^3/s)'],'.',color='mediumvioletred')
plt.xlim(melt['Date'][0],melt['Date'][2000])
plt.ylim(0,1.5)
plt.ylabel('$m^3/s$')
plt.xlabel('July 2018')
ax.xaxis.set_major_formatter(mdates.DateFormatter('%d'))
plt.savefig('dye.png',dpi=300)


fig = plt.figure(figsize=(6,2))
ax = fig.add_subplot(111)
plt.plot(dye['Time (UTC-0)'],dye['Discharge (m^3/s)'],'.',color='mediumvioletred')
plt.plot(iceray['Time (UTC-0)'],iceray['Discharge (m^3/s)'],linewidth=1,color='darkcyan')
plt.xlim(melt['Date'][0],melt['Date'][2000])
plt.ylim(0,1.5)
plt.ylabel('$m^3/s$')
plt.xlabel('July 2018')
ax.xaxis.set_major_formatter(mdates.DateFormatter('%d'))
plt.savefig('dye_iceray.png',dpi=300)

fig = plt.figure(figsize=(6,2))
ax = fig.add_subplot(111)
plt.plot(dye['Time (UTC-0)'],dye['Discharge (m^3/s)'],'.',color='mediumvioletred')
plt.plot(melt['Date'],melt['Qm3s']+0.1,linewidth=2.5,color='darkorange')
plt.plot(iceray['Time (UTC-0)'],iceray['Discharge (m^3/s)'],linewidth=1,color='darkcyan')
plt.xlim(melt['Date'][0],melt['Date'][2000])
plt.ylim(0,1.5)
plt.ylabel('$m^3/s$')
plt.xlabel('July 2018')
ax.xaxis.set_major_formatter(mdates.DateFormatter('%d'))
plt.savefig('dye_iceray_model.png',dpi=300)



