from PyMouSh import MoulinShape, TimeStamps, Qin_sinusoidal
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import pandas as pd

secinday = 24*3600
ZERO_KELVIN = 273.15

channel_length = 30000
ice_thickness = 1000

temperature_profile =np.array([ZERO_KELVIN, ZERO_KELVIN])
regional_surface_slope = 0
initial_subglacial_area = 0.1#(np.pi*0.2**2)/2)

#time
start = 0
end = 20
timestep = 300
time = np.arange(start*secinday,end*secinday,timestep)
time_figure = np.arange((start+0)*secinday,(end-0)*secinday,timestep*12)

#idealized Qin for test purposes
meltwater_input = 3
SC_lim = [0.7,0.9]
path = '/Users/celia/Dropbox/RESEARCH/MOULIN-SHAPE-FIGURES-MOVIES/MGM_movies/cylinder_5/'

# cylinder
###########################################################################
cylinder = MoulinShape(moulin_radii =5.,
                      channel_length = channel_length,
                      ice_thickness =ice_thickness)
                     
for idx,t in enumerate(time):
    cylinder.run1step(t,timestep, meltwater_input,
                    creep=False,
                    elastic_deformation=False,
                    melt_below_head=False,
                    open_channel_melt=False,
                    potential_drop=False,
                    ice_motion=False,
                    refreezing=False)
    
for idx,t_end in enumerate(time_figure):
#t_start = 10*secinday
    t_start = t_end-time_figure[-1]
    fig = plt.figure(figsize=(8,5),dpi=100)    
    cylinder.plot_MGM(fig,t_start, t_end,
                      spine_head_min=500,
                      ground_depth=-100,
                      SC_lim = SC_lim,
                      Q_fixed = True) 

    # plt.savefig(path+'cylinder5_mgm_%d.png'%idx)
    # plt.close(fig)


# ###########################################################################
# goblet = MoulinShape(z_elevations = [0,500,500,1000],
#                      moulin_radii = [1,1,5,5],
#                      ice_thickness =ice_thickness)
                     

