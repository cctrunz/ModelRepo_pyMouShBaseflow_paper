import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle

import seaborn as sns
sns.set()
sns.set_theme(style="ticks", font_scale=1.25)

secinday = 24*3600

folder ='compare_baseflow'

jeme_moulin = pd.read_csv('Field_Data/head_jeme.csv')
jeme_moulin = jeme_moulin.dropna()
h_real = jeme_moulin.head_bed.to_numpy()
t_real = jeme_moulin.soy.to_numpy()


picklefile = open(folder+'/bf0_fix', 'rb')
bf0_fix = pickle.load(picklefile)
picklefile.close()

picklefile = open(folder+'/bf3_fix', 'rb')
bf3_fix = pickle.load(picklefile)
picklefile.close()

picklefile = open(folder+'/bf1_osc', 'rb')
bf1_osc = pickle.load(picklefile)
picklefile.close()

picklefile = open(folder+'/bf1_osc_shift', 'rb')
bf1_osc_shift = pickle.load(picklefile)
picklefile.close()

picklefile = open(folder+'/bf_melt_shift', 'rb')
bf_melt_shift = pickle.load(picklefile)
picklefile.close()

picklefile = open(folder+'/bf1_fix', 'rb')
bf1_fix = pickle.load(picklefile)
picklefile.close()

# #%%
# fig, axs = plt.subplots(2, sharex=True, figsize=(8,5))

# axs[0].plot(bf1_osc_shift.t/secinday, bf1_osc_shift.dict['meltwater_input_moulin'], label='Qin', color='black')
# axs[0].set_ylim([0.1,0.5])
# ax0=axs[0].twinx()
# #ax0.plot(bf1_fix.t/secinday, bf1_fix.dict['subglacial_baseflow'], label='bf1_fix', color='khaki')
# #ax0.plot(bf1_osc.t/secinday, bf1_osc.dict['subglacial_baseflow'], label='bf1_osc',color='darkorange')
# ax0.plot(bf1_osc_shift.t/secinday, bf1_osc_shift.dict['subglacial_baseflow'], label='bf1_osc_shift', color='orangered')
# ax0.plot(bf_melt_shift.t/secinday, bf_melt_shift.dict['subglacial_baseflow'], label='bf_melt_shift', color='seagreen')


# ax0.spines['right'].set_color('orangered')
# ax0.yaxis.label.set_color('orangered')
# ax0.tick_params(axis='y', colors='orangered')

# #ax0.set_ylim([0.75,1.25])
# axs[0].legend(loc=2, prop={'size': 6})
# ax0.legend(loc=1, prop={'size': 6})
# axs[0].set_ylabel('$m^3/s$')

# axs[1].plot(t_real/secinday,h_real,label='field', color='black')
# #axs[1].plot(bf0_fix.t/secinday, bf0_fix.dict['head'], label='bf0_fix', color='darkgray')
# axs[1].plot(bf3_fix.t/secinday, bf3_fix.dict['head'], label='bf3_fix', color='dodgerblue')
# # axs[1].plot(bf1_fix.t/secinday, bf1_fix.dict['head'], label='bf1_fix', color='khaki')
# #axs[1].plot(bf1_osc.t/secinday, bf1_osc.dict['head'], label='bf1_osc', color='darkorange')
# axs[1].plot(bf1_osc_shift.t/secinday, bf1_osc_shift.dict['head'], label='bf1_osc_shift', color='orangered')
# axs[1].plot(bf_melt_shift.t/secinday, bf_melt_shift.dict['head'], label='bf_melt_shift', color='seagreen')


# axs[1].legend(loc=4, prop={'size': 6})

# axs[1].set_xlim([200,220]) #days
# axs[1].set_ylabel('Head (m)')
# axs[1].set_xlabel('day of year')

# plt.savefig('compare_baseflow')



#%% 

# CREATE COLOR PALETTE 

palette_lag = sns.color_palette('BuPu', n_colors=3) #RdBu PuOr
palette_mean = sns.color_palette('YlOrBr', n_colors=1) #RdBu PuOr
# sns.palplot(palette_lag)
# sns.palplot(palette_mean)

xlim = [210,230] #days
ylim_head = [100,500]
figure_size = (11,11)
size_legend = 12


#COMPARE LAG IN PEAK BASEFLOW
fig, ax = plt.subplots(3, figsize=figure_size)


ax[0].plot(bf1_osc_shift.t/secinday, bf1_osc_shift.dict['meltwater_input_moulin'], label='Meltwater input (melt model)', color='black')
ax[0].set_ylim([0.1,0.5])
ax[0].plot(bf1_osc.t/secinday, bf1_osc.dict['subglacial_baseflow'], label='Baseflow, 1  $m^3/s$',color=palette_lag[0])
ax[0].plot(bf1_fix.t/secinday, bf1_fix.dict['subglacial_baseflow'], label='Baseflow, 1 +/- 0.1 $m^3/s$', color=palette_lag[1])
ax[0].plot(bf1_osc_shift.t/secinday, bf1_osc_shift.dict['subglacial_baseflow'], label='Baseflow, 1 +/- 0.1 $m^3/s$, 6h lag ', color=palette_lag[2])
ax[0].plot(bf3_fix.t/secinday, bf3_fix.dict['subglacial_baseflow'], label='Baseflow, 3$m^3/s$',color=palette_mean[0])


#plot preferences
ax[0].legend(loc=1, prop={'size': size_legend})
ax[0].set_xlim(xlim) 
ax[0].set_ylim([0,3.5]) 
ax[0].set_ylabel('$m^3/s$')


#plot base head 
ax[1].plot(t_real/secinday,h_real,label='Field data', color='black')
ax[1].plot(bf0_fix.t/secinday, bf0_fix.dict['head'], label='No baseflow', color='darkgray')
#plot head with lag
ax[1].plot(bf1_osc.t/secinday, bf1_osc.dict['head'], color=palette_lag[0]) # label='oscillating baseflow in phase'
ax[1].plot(bf1_fix.t/secinday, bf1_fix.dict['head'], color=palette_lag[1]) # label='fixed baseflow'
ax[1].plot(bf1_osc_shift.t/secinday, bf1_osc_shift.dict['head'], color=palette_lag[2]) # label='oscillating baseflow out of phase'
#plot preferences
ax[1].legend(loc=3, prop={'size': size_legend})
ax[1].set_xlim(xlim) 
ax[1].set_ylim(ylim_head)
ax[1].set_ylabel('Head (m)')




#COMPARE mean BASEFLOW

#plot base head 
ax[2].plot(t_real/secinday,h_real, color='black') # label='Field data',
ax[2].plot(bf0_fix.t/secinday, bf0_fix.dict['head'], color='darkgray') # label='no baseflow', 
#plot head with lag
ax[2].plot(bf1_osc_shift.t/secinday, bf1_osc_shift.dict['head'], color=palette_lag[2]) # label='oscillating baseflow out of phase'
ax[2].plot(bf3_fix.t/secinday, bf3_fix.dict['head'], color=palette_mean[0]) # label='fixed baseflow (3$m^3/s$)'

#plot preferences
ax[2].set_xlim(xlim) 
ax[2].set_ylim(ylim_head)
ax[2].set_ylabel('Head (m)')
ax[2].set_xlabel('day of year')
xlabels = ['{:,.1f}'.format(x) for x in ax[2].get_xticks()]
ax[2].set_xticklabels(xlabels)



#COMPARE baseflow and Qin

sns.despine(offset=10,trim=True)
sns.despine(ax=ax[0],bottom=True)
sns.despine(ax=ax[1],bottom=True)
ax[0].set_xticks([])
ax[1].set_xticks([])



plt.savefig('compare_baseflow_lag.pdf')




































