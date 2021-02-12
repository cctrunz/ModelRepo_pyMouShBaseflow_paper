import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle
secinday = 24*3600

jeme_moulin = pd.read_csv('Field_Data/head_jeme.csv')
jeme_moulin = jeme_moulin.dropna()
h_real = jeme_moulin.head_bed.to_numpy()
t_real = jeme_moulin.soy.to_numpy()


picklefile = open('bf0_fix', 'rb')
bf0_fix = pickle.load(picklefile)
picklefile.close()

picklefile = open('bf3_fix', 'rb')
bf3_fix = pickle.load(picklefile)
picklefile.close()

picklefile = open('bf1_osc', 'rb')
bf1_osc = pickle.load(picklefile)
picklefile.close()

picklefile = open('bf1_osc_shift', 'rb')
bf1_osc_shift = pickle.load(picklefile)
picklefile.close()

picklefile = open('bf_melt_shift', 'rb')
bf_melt_shift = pickle.load(picklefile)
picklefile.close()

picklefile = open('bf1_fix', 'rb')
bf1_fix = pickle.load(picklefile)
picklefile.close()

#%%
fig, axs = plt.subplots(2, sharex=True, figsize=(7,4))
fig.suptitle('Vertically stacked subplots')

axs[0].plot(bf1_osc_shift.t/secinday, bf1_osc_shift.dict['meltwater_input_moulin'], label='Qin', color='black')
axs[0].set_ylim([0.1,0.5])
ax0=axs[0].twinx()
ax0.plot(bf1_fix.t/secinday, bf1_fix.dict['subglacial_baseflow'], label='bf1_fix', color='khaki')
ax0.plot(bf1_osc.t/secinday, bf1_osc.dict['subglacial_baseflow'], label='bf1_osc',color='darkorange')
ax0.plot(bf1_osc_shift.t/secinday, bf1_osc_shift.dict['subglacial_baseflow'], label='bf1_osc_shift', color='orangered')
ax0.plot(bf_melt_shift.t/secinday, bf_melt_shift.dict['subglacial_baseflow'], label='bf_melt_shift', color='seagreen')


ax0.spines['right'].set_color('orangered')
ax0.yaxis.label.set_color('orangered')
ax0.tick_params(axis='y', colors='orangered')

ax0.set_ylim([0.75,1.25])
axs[0].legend()
ax0.legend()
axs[0].set_ylabel('$m^3/s$')

axs[1].plot(t_real/secinday,h_real,label='field', color='black')
axs[1].plot(bf0_fix.t/secinday, bf0_fix.dict['head'], label='bf0_fix', color='darkgray')
axs[1].plot(bf3_fix.t/secinday, bf3_fix.dict['head'], label='bf3_fix', color='dodgerblue')
axs[1].plot(bf1_osc.t/secinday, bf1_fix.dict['head'], label='bf1_fix', color='khaki')
axs[1].plot(bf1_osc.t/secinday, bf1_osc.dict['head'], label='bf1_osc', color='darkorange')
axs[1].plot(bf1_osc_shift.t/secinday, bf1_osc_shift.dict['head'], label='bf1_osc_shift', color='orangered')
axs[1].plot(bf_melt_shift.t/secinday, bf_melt_shift.dict['head'], label='bf_melt_shift', color='seagreen')


axs[1].legend()

axs[1].set_xlim([210,220]) #days
axs[1].set_ylabel('Head (m)')
axs[1].set_xlabel('day of year')