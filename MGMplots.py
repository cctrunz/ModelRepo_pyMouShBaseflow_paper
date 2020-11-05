from PyMouSh import MoulinShape, TimeStamps, Qin_sinusoidal
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import pandas as pd
import os

secinday = 24*3600
ZERO_KELVIN = 273.15

channel_length = 30000
ice_thickness = 1000
heq = 707.38

temperature_profile =np.array([ZERO_KELVIN, ZERO_KELVIN])
regional_surface_slope = 0
initial_subglacial_area = 0.1#(np.pi*0.2**2)/2)

#time
start = 0
end = 20
timestep = 300
time = np.arange(start*secinday,end*secinday,timestep)
time_figure = np.arange((start+5)*secinday,(end-5)*secinday,timestep)

#idealized Qin for test purposes
Qin_mean = 3
dQ = 0.4 
meltwater_input = Qin_sinusoidal(time,Qin_mean, dQ)
Q_lim = [2,4]
SC_lim = [0.85,0.95]
#path = '/Users/celia/Dropbox/RESEARCH/MOULIN-SHAPE-FIGURES-MOVIES/MGM_movies/' # for little macbook pro
#path = '/Users/cctrunz/Dropbox/RESEARCH/MOULIN-SHAPE-FIGURES-MOVIES/MGM_movies/' # for big mac
path = 'C:/Users/celia/Dropbox/RESEARCH/MOULIN-SHAPE-FIGURES-MOVIES/MGM_movies/' # for surface pro

# cylinder = MoulinShape(moulin_radii = 1.,
#                         channel_length = channel_length,
#                         ice_thickness =ice_thickness)
                     
# for idx,t in enumerate(time):
#     cylinder.run1step(t,timestep, meltwater_input[idx],
#                     creep=False,
#                     elastic_deformation=False,
#                     melt_below_head=False,
#                     open_channel_melt=False,
#                     potential_drop=False,
#                     ice_motion=False,
#                     refreezing=False)
# t_start = 0
# t_end = t_start + 5*secinday
# fig = plt.figure(figsize=(6,4),dpi=300)    
# cylinder.plot_MGM(fig,t_start, t_end,
#                   spine_head_min=500,
#                   ground_depth=-100,
#                   Q_lim = Q_lim,
#                   SC_lim = SC_lim,
#                   Q_fixed = False)
    
idx_dir = 12  
   
#%%##########################################################################

#for r_top in [1.,4.,8.]:
r_top = 8
for r_bottom in [1.,4.,8.]:
    for z_break in [500,heq]:
        idx_dir = idx_dir+1
        os.mkdir('C:/Users/celia/Dropbox/RESEARCH/MOULIN-SHAPE-FIGURES-MOVIES/MGM_movies/folder_%d'%idx_dir)
        moulin = MoulinShape(z_elevations = [0,z_break,z_break,1000],
                              moulin_radii = [r_bottom,r_bottom,r_top,r_top],
                              channel_length = channel_length,
                              ice_thickness = ice_thickness)
                             
        for idx,t in enumerate(time):
            moulin.run1step(t,timestep, meltwater_input[idx],
                            creep=False,
                            elastic_deformation=False,
                            melt_below_head=False,
                            open_channel_melt=False,
                            potential_drop=False,
                            ice_motion=False,
                            refreezing=False)
         
        for idx,t_start in enumerate(time_figure):
            t_end = t_start + 5*secinday
            fig = plt.figure(figsize=(6,4))    
            moulin.plot_MGM(fig,t_start, t_end,
                              spine_head_min=500,
                              ground_depth=-100,
                              Q_lim = Q_lim,
                              SC_lim = SC_lim,
                              Q_fixed = False)
            plt.savefig(path+'folder_%d/MGM_top%d_bottom%d_zbreak%d_no%d.png'%(idx_dir,r_top,r_bottom,z_break,idx))#,dpi=200)
            plt.clf()
            plt.close(fig)
            
    del moulin
            
#%%##########################################################################
for middle in [1.,4.,8.]:
    for main in [1.,4.,8.]:
        idx_dir = idx_dir+1
        os.mkdir('C:/Users/celia/Dropbox/RESEARCH/MOULIN-SHAPE-FIGURES-MOVIES/MGM_movies/folder_%d'%idx_dir)
        losange = MoulinShape(z_elevations = [0,heq-(1000-heq),heq,1000],
                              moulin_radii = [main,main,middle,main],
                              channel_length = channel_length,
                              ice_thickness =ice_thickness)
                             
        for idx,t in enumerate(time):
            losange.run1step(t,timestep, meltwater_input[idx],
                            creep=False,
                            elastic_deformation=False,
                            melt_below_head=False,
                            open_channel_melt=False,
                            potential_drop=False,
                            ice_motion=False,
                            refreezing=False)
         
        for idx,t_start in enumerate(time_figure):
            t_end = t_start + 5*secinday
            fig = plt.figure(figsize=(8,6))   
            losange.plot_MGM(fig,t_start, t_end,
                              spine_head_min=500,
                              ground_depth=-100,
                              Q_lim = Q_lim,
                              SC_lim = SC_lim,
                              Q_fixed = False)
            plt.savefig(path+'folder_%dLosange_middle%d_main%d_no%d.png'%(idx_dir,middle,main,idx))
            plt.clf()
            plt.close(fig)    
    del losange
