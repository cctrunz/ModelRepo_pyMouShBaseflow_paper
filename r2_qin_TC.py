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
dye = dye_raw.Q.groupby(pd.Grouper(freq='1H')).mean()
#stage = stage_raw.stage.groupby(pd.Grouper(freq='1H')).mean()

#
#stage_long = stage[(stage.index > stage.index[0]) & (stage.index < stage.index[-1])]  
#model_long = model_raw.Q[(model_raw.index > dye.index[0]) & (model_raw.index < dye.index[-1])]       

# same period of time
#stage = stage[(stage.index > dye.index[0]) & (stage.index < dye.index[-1])]  
model = model_raw.Q[(model_raw.index > dye.index[0]) & (model_raw.index < dye.index[-1])]     

#add missing timesteps and fill missing days without data with nans to obtain homogenous data
#idx_long = pd.date_range(stage_long.index[0], stage_long.index[-1], freq='1H')
#stage_long = stage_long.reindex(idx_long, fill_value=np.NaN)
#model_long = model_long.reindex(idx_long, fill_value=np.NaN)

#add missing timesteps and fill missing days without data with nans to obtain homogenous data
idx = pd.date_range(dye.index[0], dye.index[-1], freq='1H')
#stage = stage.reindex(idx, fill_value=np.NaN)
dye = dye.reindex(idx, fill_value=np.NaN)
model = model.reindex(idx, fill_value=np.NaN)

#%% Comparison with dye
date_mid = '2017-07-24 07:00:00'
# varx = dye.values
# varx_a = dye[dye.index[0]:date_mid].values
# varx_b = dye[date_mid:dye.index[-1]].values
# varys =  stage.values
# varys_a =  stage[dye.index[0]:date_mid].values
# varys_b =  stage[date_mid:dye.index[-1]].values
# varym =  model.values
# varym_a =  model[dye.index[0]:date_mid].values
# varym_b =  model[date_mid:dye.index[-1]].values

# varx = stage.values  
# varx_a = stage[dye.index[0]:date_mid].values
# varx_b = stage[date_mid:dye.index[-1]].values

varx =  dye[1:-2].values
# varys_a =   dye[dye.index[0]:date_mid].values
# varys_b =   dye[date_mid:dye.index[-1]].values

vary =  model[1:-2].values
# varym_a =  model[dye.index[0]:date_mid].values
# varym_b =  model[date_mid:dye.index[-1]].values

# mask_stage = ~np.isnan(varx) & ~np.isnan(varys)
# mask_stage_a = ~np.isnan(varx_a) & ~np.isnan(varys_a)
# mask_stage_b = ~np.isnan(varx_b) & ~np.isnan(varys_b)

# mask_model = ~np.isnan(varx) & ~np.isnan(varym)
# mask_model_a = ~np.isnan(varx_a) & ~np.isnan(varym_a)
# mask_model_b = ~np.isnan(varx_b) & ~np.isnan(varym_b)

# slope_stage, intercept_stage, r_value_stage, p_value_stage, std_err_stage = stats.linregress(varx[mask_stage], varys[mask_stage])
# slope_stage_a, intercept_stage_a, r_value_stage_a, p_value_stage_a, std_err_stage_a = stats.linregress(varx_a[mask_stage_a], varys_a[mask_stage_a])
# slope_stage_b, intercept_stage_b, r_value_stage_b, p_value_stage_b, std_err_stage_b = stats.linregress(varx_b[mask_stage_b], varys_b[mask_stage_b])

slope, intercept, r_value, p_value, std_err = stats.linregress(varx, vary)
# slope_model_a, intercept_model_a, r_value_model_a, p_value_model_a, std_err_model_a = stats.linregress(varx_a[mask_model_a], varym_a[mask_model_a])
# slope_model_b, intercept_model_b, r_value_model_b, p_value_model_b, std_err_model_b = stats.linregress(varx_b[mask_model_b], varym_b[mask_model_b])

import math



RMSE = math.sqrt(np.square(np.subtract(vary,varx)).mean())

print("Root Mean Square Error:\n")
print(RMSE)


#%%

n = len(dye)
color=cm.rainbow(np.linspace(0,1,n))

plt.figure()
for i,c in zip(range(n),color):
    plt.plot(dye.values[i],model.values[i],marker='o', linestyle='', c=c)
    plt.plot(varx,varx*slope+intercept,color='grey')

#plt.plot(dye.values,stage.values,marker='o', linestyle='')
plt.text(2.85,0.3,'All (grey): R$^2$=%s, p=%s\nDay1'
         %(round(r_value**2,2),p_value
           ))
plt.xlabel('Discharge (dye)')
plt.ylabel('Discharge (model)')



#%%
#%%
#%%

# all on the same timestep
dye = dye_raw.Q.groupby(pd.Grouper(freq='1H')).mean()
stage = stage_raw.stage.groupby(pd.Grouper(freq='1H')).mean()

#
stage_long = stage[(stage.index > stage.index[0]) & (stage.index < stage.index[-1])]  
model_long = model_raw.Q[(model_raw.index > stage.index[0]) & (model_raw.index < stage.index[-1])]       

# same period of time
stage = stage[(stage.index > dye.index[0]) & (stage.index < dye.index[-1])]  
model = model_raw.Q[(model_raw.index > dye.index[0]) & (model_raw.index < dye.index[-1])]     

#add missing timesteps and fill missing days without data with nans to obtain homogenous data
idx_long = pd.date_range(stage_long.index[0], stage_long.index[-1], freq='1H')
stage_long = stage_long.reindex(idx_long, fill_value=np.NaN)
model_long = model_long.reindex(idx_long, fill_value=np.NaN)

#add missing timesteps and fill missing days without data with nans to obtain homogenous data
idx = pd.date_range(dye.index[0], dye.index[-1], freq='1H')
stage = stage.reindex(idx, fill_value=np.NaN)
dye = dye.reindex(idx, fill_value=np.NaN)
model = model.reindex(idx, fill_value=np.NaN)

#%% Comparison with dye
date_mid = '2017-07-24 07:00:00'
# varx = dye.values
# varx_a = dye[dye.index[0]:date_mid].values
# varx_b = dye[date_mid:dye.index[-1]].values
# varys =  stage.values
# varys_a =  stage[dye.index[0]:date_mid].values
# varys_b =  stage[date_mid:dye.index[-1]].values
# varym =  model.values
# varym_a =  model[dye.index[0]:date_mid].values
# varym_b =  model[date_mid:dye.index[-1]].values

varx = stage.values  
varx_a = stage[dye.index[0]:date_mid].values
varx_b = stage[date_mid:dye.index[-1]].values

varys =  dye.values
varys_a =   dye[dye.index[0]:date_mid].values
varys_b =   dye[date_mid:dye.index[-1]].values

varym =  model.values
varym_a =  model[dye.index[0]:date_mid].values
varym_b =  model[date_mid:dye.index[-1]].values

mask_stage = ~np.isnan(varx) & ~np.isnan(varys)
mask_stage_a = ~np.isnan(varx_a) & ~np.isnan(varys_a)
mask_stage_b = ~np.isnan(varx_b) & ~np.isnan(varys_b)

mask_model = ~np.isnan(varx) & ~np.isnan(varym)
mask_model_a = ~np.isnan(varx_a) & ~np.isnan(varym_a)
mask_model_b = ~np.isnan(varx_b) & ~np.isnan(varym_b)

slope_stage, intercept_stage, r_value_stage, p_value_stage, std_err_stage = stats.linregress(varx[mask_stage], varys[mask_stage])
slope_stage_a, intercept_stage_a, r_value_stage_a, p_value_stage_a, std_err_stage_a = stats.linregress(varx_a[mask_stage_a], varys_a[mask_stage_a])
slope_stage_b, intercept_stage_b, r_value_stage_b, p_value_stage_b, std_err_stage_b = stats.linregress(varx_b[mask_stage_b], varys_b[mask_stage_b])

slope_model, intercept_model, r_value_model, p_value_model, std_err_model = stats.linregress(varx[mask_model], varym[mask_model])
slope_model_a, intercept_model_a, r_value_model_a, p_value_model_a, std_err_model_a = stats.linregress(varx_a[mask_model_a], varym_a[mask_model_a])
slope_model_b, intercept_model_b, r_value_model_b, p_value_model_b, std_err_model_b = stats.linregress(varx_b[mask_model_b], varym_b[mask_model_b])

import math



RMSE_model = math.sqrt(np.square(np.subtract(varym[mask_model],varx[mask_model])).mean())
RMSE_model_a = math.sqrt(np.square(np.subtract(varym_a[mask_model_a],varx_a[mask_model_a])).mean())
RMSE_model_b = math.sqrt(np.square(np.subtract(varym_b[mask_model_b],varx_b[mask_model_b])).mean())
print("Root Mean Square Error:\n")
print(RMSE_model, RMSE_model_a, RMSE_model_b)




n = len(dye)
color=cm.rainbow(np.linspace(0,1,n))

plt.figure()
for i,c in zip(range(n),color):
    plt.plot(stage.values[i],dye.values[i],marker='o', linestyle='', c=c)
    plt.plot(varx,varx*slope_stage+intercept_stage,color='grey')
    plt.plot(varx_a,varx_a*slope_stage_a+intercept_stage_a,color=color[round(n/4)])
    plt.plot(varx_b,varx_b*slope_stage_b+intercept_stage_b,color=color[round(3*n/4)])
#plt.plot(dye.values,stage.values,marker='o', linestyle='')
plt.text(2.85,0.3,'All (grey): R$^2$=%s, p=%s\nDay1 (blue): R$^2$=%s, p=%s\nDay2(orange): R$^2$=%s, p=%s'
         %(round(r_value_stage**2,2),p_value_stage,
           round(r_value_stage_a**2,2),p_value_stage_a,
           round(r_value_stage_b**2,2),p_value_stage_b
           ))
plt.xlabel('Stream water level')
plt.ylabel('Discharge (dye)')

plt.figure()
for i,c in zip(range(n),color):
    plt.plot(stage.values[i],model.values[i],marker='o', linestyle='', c=c)
    plt.plot(varx,varx*slope_model+intercept_model,color='grey')
    plt.plot(varx_a,varx_a*slope_model_a+intercept_model_a,color=color[round(n/4)])
    plt.plot(varx_b,varx_b*slope_model_b+intercept_model_b,color=color[round(3*n/4)])
#plt.plot([0,0.35],[0,0.35], linestyle='--')
plt.text(2.85,0.2,'All (grey): R$^2$=%s, p=%s\nDay1 (blue): R$^2$=%s, p=%s\nDay2(orange): R$^2$=%s, p=%s'
         %(round(r_value_model**2,2),p_value_model,
           round(r_value_model_a**2,2),p_value_model_a,
           round(r_value_model_b**2,2),p_value_model_b
           ))
plt.xlabel('Stream water level')
plt.ylabel('Discharge (model)')



# #%% Comparison dye vs model

# varx = dye.values

# varys =  stage.values
# varys_a =  stage[dye.index[0]:date_mid].values
# varys_b =  stage[date_mid:dye.index[-1]].values
# varym =  model.values
# varym_a =  model[dye.index[0]:date_mid].values
# varym_b =  model[date_mid:dye.index[-1]].values
 
# mask_stage = ~np.isnan(varx) & ~np.isnan(varys)
# mask_stage_a = ~np.isnan(varx_a) & ~np.isnan(varys_a)
# mask_stage_b = ~np.isnan(varx_b) & ~np.isnan(varys_b)

# mask_model = ~np.isnan(varx) & ~np.isnan(varym)
# mask_model_a = ~np.isnan(varx_a) & ~np.isnan(varym_a)
# mask_model_b = ~np.isnan(varx_b) & ~np.isnan(varym_b)

# slope_stage, intercept_stage, r_value_stage, p_value_stage, std_err_stage = stats.linregress(varx[mask_stage], varys[mask_stage])
# slope_stage_a, intercept_stage_a, r_value_stage_a, p_value_stage_a, std_err_stage_a = stats.linregress(varx_a[mask_stage_a], varys_a[mask_stage_a])
# slope_stage_b, intercept_stage_b, r_value_stage_b, p_value_stage_b, std_err_stage_b = stats.linregress(varx_b[mask_stage_b], varys_b[mask_stage_b])

# slope_model, intercept_model, r_value_model, p_value_model, std_err_model = stats.linregress(varx[mask_model], varym[mask_model])
# slope_model_a, intercept_model_a, r_value_model_a, p_value_model_a, std_err_model_a = stats.linregress(varx_a[mask_model_a], varym_a[mask_model_a])
# slope_model_b, intercept_model_b, r_value_model_b, p_value_model_b, std_err_model_b = stats.linregress(varx_b[mask_model_b], varym_b[mask_model_b])

# import math



# RMSE_model = math.sqrt(np.square(np.subtract(varym[mask_model],varx[mask_model])).mean())
# RMSE_model_a = math.sqrt(np.square(np.subtract(varym_a[mask_model_a],varx_a[mask_model_a])).mean())
# RMSE_model_b = math.sqrt(np.square(np.subtract(varym_b[mask_model_b],varx_b[mask_model_b])).mean())
# print("Root Mean Square Error:\n")
# print(RMSE_model, RMSE_model_a, RMSE_model_b)

#%%
# plt.figure()
# plt.plot(model)
# plt.plot(stage)
# plt.plot(dye)

# n = len(dye)
# color=cm.rainbow(np.linspace(0,1,n))

# plt.figure()
# for i,c in zip(range(n),color):
#     plt.plot(dye.values[i],stage.values[i],marker='o', linestyle='', c=c)
#     plt.plot(varx,varx*slope_stage+intercept_stage,color='grey')
#     plt.plot(varx_a,varx_a*slope_stage_a+intercept_stage_a,color=color[round(n/4)])
#     plt.plot(varx_b,varx_b*slope_stage_b+intercept_stage_b,color=color[round(3*n/4)])
# #plt.plot(dye.values,stage.values,marker='o', linestyle='')
# plt.text(0.05,3.16,'All (grey): R$^2$=%s, p=%s\nDay1 (blue): R$^2$=%s, p=%s\nDay2(orange): R$^2$=%s, p=%s'
#          %(round(r_value_stage**2,2),p_value_stage,
#            round(r_value_stage_a**2,2),p_value_stage_a,
#            round(r_value_stage_b**2,2),p_value_stage_b
#            ))
# plt.xlabel('Discharge (dye)')
# plt.ylabel('Stream water level')

# plt.figure()
# for i,c in zip(range(n),color):
#     plt.plot(dye.values[i],model.values[i],marker='o', linestyle='', c=c)
#     plt.plot(varx,varx*slope_model+intercept_model,color='grey')
#     plt.plot(varx_a,varx_a*slope_model_a+intercept_model_a,color=color[round(n/4)])
#     plt.plot(varx_b,varx_b*slope_model_b+intercept_model_b,color=color[round(3*n/4)])
# plt.plot([0,0.35],[0,0.35], linestyle='--')
# plt.text(0.,0.29,'All (grey): R$^2$=%s, p=%s\nDay1 (blue): R$^2$=%s, p=%s\nDay2(orange): R$^2$=%s, p=%s'
#          %(round(r_value_model**2,2),p_value_model,
#            round(r_value_model_a**2,2),p_value_model_a,
#            round(r_value_model_b**2,2),p_value_model_b
#            ))
# plt.xlabel('Discharge (dye)')
# plt.ylabel('Discharge (model)')

#%%
# plt.figure()
# for i,c in zip(range(n),color):
#     plt.plot(stage_long.values[i],model_long.values[i],marker='o', linestyle='', c=c)
#     #plt.plot(varys_long,varys_long*slope_long+intercept_long,color='grey')
# #plt.plot(dye.values,stage.values,marker='o', linestyle='')
# plt.text(0.05,3.16,'All (grey): R$^2$=%s, p=%s'
#          %(round(r_value_long**2,2),p_value_long
#            ))
# plt.xlabel('Discharge (dye)')
# plt.ylabel('Stream water level')



#