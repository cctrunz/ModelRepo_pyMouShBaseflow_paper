#import datetime
#from datetime import datetime as dt 

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#import matplotlib.dates as mdates
#import matplotlib.units as munits

# format matplotlib dates
# converter = mdates.ConciseDateConverter()
# munits.registry[np.datetime64] = converter
# munits.registry[datetime.date] = converter
# munits.registry[dt] = converter


# define functions to implement SUH equations
def synthetic_unit_hydrograph(m, tp, hp):
    """create the SUH q-curve with empirical values for each catchment. 
    
    SUH: synthetic unit hydrograph
        update SUH with empirical values m, tp, and hp
    Implements the equation used by King (2018) pp 41.
    
    Parameters
    ----------
        m (float, np.array) : empirical equation shape factor
        tp (np.array) : time to peak discharge (hrs)
        hp (np.array): peak discharge (hr^-1)
        
    Returns
    -------
        q (np.array) : q-curve units of m/hr
    """
    return (np.exp(m) * (np.arange(25) / tp) ** m 
            * (np.exp( -m * (np.arange(25) / tp))) * hp )

def calc_QSUH(melt, db=None, c=None, m=None, tp=None, hp=None):
    """Calculate discharge using the synthetic unit hydrograph for melt timeseries.
    
    Parameters
    ----------
        melt (pd.Series) : melt timeseries
        const (np.array, pd.df, pd.series) : constants needed for SUH
    
    
    Returns
    -------
        Q (np.array) : Discharge calculated with SUH (m^3/s)
        
    """
    #param = param_SUH.loc[moulin_name]

    #calc UH
    # M' effective 
    # units: mm/h
    runoff = melt* c 
    
    # Synthetic Unit Hydrograph - King (2018)
    SUH = synthetic_unit_hydrograph(m, 
                                    tp, 
                                    hp)
    # convolution of runoff and the unit hydrograph q-curve
    # Q = M' * q where * is the convolution operator (np.convolve)
    q = np.convolve(runoff,SUH)
    # remove the adds 23 additional values
    q = q[0:-24] 
    # multiply by drainage basin area and convert 
    # runoff from mm per hour to cubic meters per second
    # <Note!> Before this point, units are in 1/hours.
    Q = q * db*1e3/3600
    Q = pd.Series(data=Q,index=runoff.index).rename('Discharge (m^3/s)')
    
    # if visualize:
    #     gen_figs(SUH, Q, moulin_name)
    
    return SUH,Q


#import low camp melt and weather
lowc = pd.read_csv('LOWC_MELT_2017_2018.csv', 
                    index_col=0, 
                    parse_dates=True)

# import discharge measurements for plotting and SUH calibration
Q_JEME = pd.read_csv('JEME_Qdye.csv', index_col=0, parse_dates=True)
Q_PIRA = pd.read_csv('PIRA_Q.csv', index_col=0, parse_dates=True, header=[0,1])


param_SUH = pd.read_csv('param_SUH.csv',index_col=0)
param_SUH['hp'] = np.round(param_SUH['Cp']/param_SUH['tp'],3)
# param_SUH


#%%
#import SUH parameters
# param_SUH = {'db': [0.24],#[0.2,0.5], 
#              'c':  [0.4,0.9],
#              'Cp': [0.6,1],
#              'm':  [2.5,1],
#              'tp': [2.6]}




param_jeme = param_SUH.loc['JEME']
db_list = [param_SUH.loc['JEME'].db-param_SUH.loc['JEME'].db*0.5, 
           param_SUH.loc['JEME'].db,
           param_SUH.loc['JEME'].db+param_SUH.loc['JEME'].db*0.5]

c_list = [0.5, 0.9]
m_list = [1, 2, 4]
Cp_list = [0.5, 0.6, 0.7]
tp_value= [param_jeme.tp]



    
listdic =[]

for db, c, Cp, m, tp in [(db,c,Cp,m,tp)  for db in db_list
                                          for c in c_list
                                          for Cp in Cp_list
                                          for m in m_list
                                          for tp in tp_value]:
    
    
    hp = np.round(Cp/tp,3)
    print('db=', db,' / c=', c, ' / Cp=', Cp, ' / m=', m, ' / tp=', tp) 
    
    SUH, Q = calc_QSUH(lowc.melt_rate, db=db, c=c, m=m, tp=tp, hp=hp)
    
    results = {'db': db, 
                'c': c,
                'Cp': Cp,
                'm': m,
                'tp': tp,
                'hp': hp,
                'SUH': SUH,
                'Q': Q}
    
    listdic.append(results)
    
#import SUH parameters
#param_SUH = pd.read_csv('param_SUH.csv',index_col=0)
# param_SUH.at['JEME','Cp'] = 1
# param_SUH.at['PIRA','Cp'] = 1
# param_SUH.at['JEME','db'] = 0.3
# param_SUH.at['PIRA','db'] = 0.5
# param_SUH.at['JEME','c'] = 0.9
# param_SUH.at['PIRA','c'] = 0.9
# param_SUH.at['JEME','m'] = 1
# param_SUH.at['PIRA','m'] = 1
#%%

#I should do my choices, +/- the range of possibilites

#compare drainage basins db
#############################

fig=plt.figure(tight_layout=True, figsize=(5,6))
#ax1 = fig.add_subplot(311, ylabel='q', xlabel='hours')
ax1 = fig.add_subplot(411, ylabel='Q m$^3$/s')
ax2 = fig.add_subplot(412, ylabel='Q m$^3$/s')
ax3 = fig.add_subplot(413, ylabel='Q m$^3$/s')
ax4 = fig.add_subplot(414, ylabel='Q m$^3$/s')


list_select = list(filter(lambda dic: dic['db'] == param_jeme['db']   
                            and dic['c'] == param_jeme['c']
                            and dic['Cp'] == param_jeme['Cp']
                            and dic['m'] == param_jeme['m'], 
                            listdic))

list_db = list(filter(lambda dic: dic['c'] == param_jeme['c']
                            and dic['Cp'] == param_jeme['Cp']
                            and dic['m'] == param_jeme['m'], 
                            listdic))

list_c = list(filter(lambda dic: dic['db'] == param_jeme['db']                    
                            and dic['Cp'] == param_jeme['Cp']
                            and dic['m'] == param_jeme['m'], 
                            listdic))

list_Cp = list(filter(lambda dic: dic['db'] == param_jeme['db']
                            and dic['c'] == param_jeme['c']
                            and dic['m'] == param_jeme['m'], 
                            listdic))

list_m = list(filter(lambda dic: dic['db'] == param_jeme['db']
                            and dic['c'] == param_jeme['c']
                            and dic['Cp'] == param_jeme['Cp'], 
                            listdic))
#ax1.title(str(list_db))
#♫ax1.text(15,0.1,'DB = %s \n c = %s\n m = %s \n Cp = %s \n tp = %s'%(list_filtered[0]['db'], list_filtered[0]['c'], list_filtered[0]['m'], list_filtered[0]['Cp'], list_filtered[0]['tp']))

# for dic in list_db:   
#     #if dic['db']
#     ax1.plot(np.arange(25),dic['SUH'])
#     ax2.plot(dic['Q'],color='grey')
#     ax3.plot(dic['Q'],color='grey')

#Plot Db variation
color='grey'
alpha=0.3
ax1.fill_between(list_db[0]['Q'].index, list_db[0]['Q'], list_db[-1]['Q'], color=color, alpha=alpha)
ax2.fill_between(list_c[0]['Q'].index, list_c[0]['Q'], list_c[1]['Q'], color=color, alpha=alpha)
ax3.fill_between(list_Cp[0]['Q'].index, list_Cp[0]['Q'], list_Cp[1]['Q'], color=color, alpha=alpha)
ax4.fill_between(list_m[0]['Q'].index, list_m[0]['Q'], list_m[1]['Q'], color=color, alpha=alpha)



#Plot real values    


ax1.plot(Q_JEME['Discharge_m^3/s'],'.',markersize=1, color='blue', label='dye, measured')
ax1.set_xlim(Q_JEME.index[0]-pd.Timedelta('12H'),Q_JEME.index[-1]+pd.Timedelta('12H'))
ax1.plot(list_select[0]['Q'], color='black')

ax2.plot(Q_JEME['Discharge_m^3/s'],'.',markersize=1, color='blue', label='dye, measured')
ax2.set_xlim(Q_JEME.index[0]-pd.Timedelta('12H'),Q_JEME.index[-1]+pd.Timedelta('12H'))
ax2.plot(list_select[0]['Q'], color='black')

ax3.plot(Q_JEME['Discharge_m^3/s'],'.',markersize=1, color='blue', label='dye, measured')
ax3.set_xlim(Q_JEME.index[0]-pd.Timedelta('12H'),Q_JEME.index[-1]+pd.Timedelta('12H'))
ax3.plot(list_select[0]['Q'], color='black')

ax4.plot(Q_JEME['Discharge_m^3/s'],'.',markersize=1, color='blue', label='dye, measured')
ax4.set_xlim(Q_JEME.index[0]-pd.Timedelta('12H'),Q_JEME.index[-1]+pd.Timedelta('12H'))
ax4.plot(list_select[0]['Q'], color='black')


ax1.legend()

#%%

list_high = list(filter(lambda dic: dic['db'] == db_list[-1]  
                            and dic['c'] == c_list[-1]
                            and dic['Cp'] == Cp_list[-1]
                            and dic['m'] == m_list[0], 
                            listdic))

list_low = list(filter(lambda dic: dic['db'] == db_list[0]  
                            and dic['c'] == c_list[0]
                            and dic['Cp'] == Cp_list[0]
                            and dic['m'] == m_list[-1], 
                            listdic))

fig,ax = plt.subplots()

ax.fill_between(list_high[0]['Q'].index,list_low[0]['Q'],list_high[0]['Q'], color='lightgrey', label='uncertainty')
ax.plot(Q_JEME['Discharge_m^3/s'],'.',markersize=1, color='blue', label='dye, measured')
ax.plot(Q_PIRA['ice_ray','Discharge (m^3/s)'],'.',markersize=1, color='blue')
#ax.plot(Q, label='modeled discharge')
ax.plot(list_select[0]['Q'], color='black', label='modeled')
ax.set_xlim(Q_JEME.index[0]-pd.Timedelta('12H'),Q_JEME.index[-1]+pd.Timedelta('12H'))
ax.set_ylim(0,0.5)
ax.legend()

#ax.set_xlim(Q_PIRA.index[0]-pd.Timedelta('12H'),Q_PIRA.index[-1]+pd.Timedelta('12H'))
plt.savefig('blind_estimate_jeme.png',dpi=300) 


#%%

runoff = 0.9

Q_low = lowc.melt_rate *db_list[0]*1e3/3600 * runoff
Q = lowc.melt_rate *db_list[1]*1e3/3600 * runoff
Q_high = lowc.melt_rate *db_list[2]*1e3/3600 * runoff

fig,ax = plt.subplots()

ax.fill_between(Q.index, Q_low.values, Q_high.values)
ax.plot(Q_JEME['Discharge_m^3/s'],'.',markersize=1, color='blue', label='dye, measured')
ax.plot(Q_PIRA['ice_ray','Discharge (m^3/s)'],'.',markersize=1, color='blue', label='dye, measured')
ax.plot(Q)
ax.plot(list_select[0]['Q'], color='black')
#☻ax.set_xlim(Q_JEME.index[0]-pd.Timedelta('12H'),Q_JEME.index[-1]+pd.Timedelta('12H'))
    
    
    
    
    
    