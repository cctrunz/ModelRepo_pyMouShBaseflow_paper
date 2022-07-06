#%pylab inline
from pyMouSh import MoulinShape, TimeStamps, Qin_constant, Qin_sinusoidal, Qin_real, calculate_h_S_schoof, find_nearest
import numpy as np
from datetime import datetime as dt  
import matplotlib.pylab as plt
from matplotlib.pyplot import cm
import pandas as pd
from scipy import stats
import math


### 
secinday = 24*3600

#import JEME dye
dye_raw = pd.read_csv('Field_Data/JEME_Qdye.csv', index_col=0, parse_dates=True)
#import JEME stage
stage_raw = pd.read_csv('Field_Data/LOWC17_STR_1.csv', index_col=0, parse_dates=True)

#Import meltwater input calculated with meltmodel for JEME
model_raw = pd.read_csv('Field_Data/JEME_QSUH.csv', index_col=0, parse_dates=True)


#%%





# all on the same timestep
stage = stage_raw.stage.groupby(pd.Grouper(freq='1H')).mean()


idx_long = pd.date_range(stage.index[0], stage.index[-1], freq='1H')
stage = stage.reindex(idx_long, fill_value=np.NaN)
model_raw = model_raw.reindex(idx_long, fill_value=np.NaN)
#
day_start =  '2017-07-22 16:00:00'
day_end =  '2017-07-29 16:00:00'
stage_long = stage[(stage.index > day_start) & (stage.index < day_end)]  

day_start =  '2017-07-22 19:00:00'
day_end =  '2017-07-29 19:00:00'
model_long = model_raw.Q[(model_raw.index > day_start) & (model_raw.index < day_end)]       



# Comparison with dye

varys =  stage.values

varym_long = model_long.values
varys_long = stage_long.values
mask_long  = ~np.isnan(varym_long) & ~np.isnan(varys_long)

slope_long, intercept_long, r_value_long, p_value_long, std_err_long = stats.linregress(varys_long[mask_long], varym_long[mask_long])



#
plt.figure()
plt.plot(model_long)
plt.plot(stage)


n=len(stage_long)

color=cm.rainbow(np.linspace(0,1,n))


plt.figure()
for i,c in zip(range(n),color):
    plt.plot(stage_long.values[i],model_long.values[i],marker='o', linestyle='', c=c)
    plt.plot(varys_long,varys_long*slope_long+intercept_long,color='grey')
#plt.plot(dye.values,stage.values,marker='o', linestyle='')
plt.text(2,0.3,'All (grey): R$^2$=%s, p=%s'
         %(round(r_value_long**2,2),p_value_long
           ))
plt.ylabel('Discharge (dye)')
plt.xlabel('Stream water level')



#