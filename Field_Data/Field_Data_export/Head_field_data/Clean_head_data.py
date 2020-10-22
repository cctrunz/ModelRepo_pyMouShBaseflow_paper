m3_Hbase_masl = 657 - 564
m4_Hbase_masl = 709 - 540
mf_Hbase_masl = 703 - 620



import matplotlib.pyplot as plt
from pandas import read_csv
import pandas as pd 
from datetime import datetime as dt

  

df1 = read_csv('andrews_data_1.csv')
df2 = read_csv('andrews_data_2.csv')
dfj = read_csv('JEME17_MOU_1.csv')
dfp = read_csv('PIRA18_HYD_1.csv')
dfr = read_csv('RADI17_MOU_1.csv')
dfj_idx = read_csv('JEME17_MOU_1.csv',parse_dates=True, index_col='TIMESTAMP')
dfp_idx = read_csv('PIRA18_HYD_1.csv',parse_dates=True, index_col='TIMESTAMP')
dfr_idx = read_csv('RADI17_MOU_1.csv',parse_dates=True, index_col='TIMESTAMP')



'''moulin 3'''
moulin3=pd.DataFrame()
moulin3['doy']=df1['Moulin 3 Dday']
moulin3['soy']=moulin3.doy*24*3600
moulin3['head_masl']=df1['Moulin 3 (m)']
moulin3['head_bed']=moulin3.head_masl - m3_Hbase_masl
moulin3 = moulin3.dropna()
moulin3.to_csv('head_m3.csv')

plt.figure()
plt.plot(moulin3.head_masl)

'''moulin 4'''
moulin4=pd.DataFrame()
moulin4['doy']=df1['Moulin 4 Dday']
moulin4['soy']=moulin4.doy*24*3600
moulin4['head_masl']=df1['Moulin 4 (m)']
moulin4['head_bed']=moulin4.head_masl - m3_Hbase_masl
moulin4 = moulin4.dropna()
moulin4.to_csv('head_m4.csv')

plt.figure()
plt.plot(moulin4.head_masl)

'''moulin f'''
moulinf=pd.DataFrame()
moulinf['doy']=df2['FOXX  Moulin Dday']
moulinf['soy']=moulinf.doy*24*3600
moulinf['head_masl']=df2['FOXX moulin (m)']
moulinf['head_bed']=moulinf.head_masl - m3_Hbase_masl
moulinf = moulinf.dropna()
moulinf.to_csv('head_mf.csv')

plt.figure()
plt.plot(moulinf.head_masl)


'''moulin jeme'''
jeme=pd.DataFrame()
elapsed = dfj_idx.index - dt(2017, 1, 1,00,00,00)
jeme['soy'] = elapsed.total_seconds()
jeme['doy'] = jeme.soy *24*3600
jeme['head_bed']=dfj[['water_level_above_bed']]
jeme = jeme.dropna()
jeme.to_csv('head_jeme.csv')

plt.figure()
plt.plot(jeme.head_bed)

'''moulin pira'''
pira=pd.DataFrame()
elapsed = dfp_idx.index - dt(2018, 1, 1,00,00,00)
pira['soy'] = elapsed.total_seconds()
pira['doy'] = pira.soy *24*3600
pira['head_bed']=dfp[['water_level_above_bed']]
pira = pira.dropna()
pira.to_csv('head_pira.csv')

plt.figure()
plt.plot(pira.head_bed)

#%
'''moulin radi'''
radi=pd.DataFrame()
elapsed = dfr_idx.index - dt(2017, 1, 1,00,00,00)
radi['soy'] = elapsed.total_seconds()
radi['doy'] = radi.soy *24*3600
radi['head_bed']=dfr[['water_level_above_bed']]
radi = radi.dropna()
radi.to_csv('head_radi.csv')

plt.figure()
plt.plot(radi.head_bed)

