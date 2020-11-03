''' 
Moulin Physical Model from LC Andrews and K Poinar
Translated by Celia Trunz

Component provenance:
    - Subglacial channel: Schoof2010 and Covington2020
    - Creep
    - Elastic
    - Melt
    - Refreezing
'''

import numpy as np

#from typing import Union, Tuple#, ArrayLike
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.integrate import cumtrapz
from collections import defaultdict
#import pandas as pd
#import plot_codes.comprehensive_plot

SECINDAY = 24*3600

ZERO_KELVIN = 273.15
PI = np.pi
ICE_DENSITY = 910  # kg/m3; Ice density
WATER_DENSITY = 1000  # kg/m3; Water density
GRAVITY = 9.8  # m/s2; Gravity
# J / (kg * K)   heat capacity of water for unit mass,  from Jarosch & Gundmundsson (2012)
WATER_HEAT_CAPACITY = 4210
LATENT_HEAT_FUSION = 3.32e5#335000  # J/kg
ICE_POISSON_RATIO = 0.3    # []
YOUNG_ELASTIC_MODULUS = 5e9  # Pa (Vaughan 1995) -->shearModulus in matlab
IDEAL_GAZ_CONSTANT = 8.314  # R
TEMPERATURE_TRANSITION = 263  # 10°C in kelvin
ARRHENIUS_TRANSITION = 3.5e-25  # Astar
LOW_CREEP_ACTIVATION_ENERGY = 6e4  # Qless
EFFECTIVE_CREEP_ACTIVATION_ENERGY = 11.5e4  # Qmore
ICE_EXPONENT = 3
CP = 2115  # J/kgK
KI = 2.1  # J/mKs
ICE_TEMPERATURE_DIFFUSION = KI/ICE_DENSITY/CP
# A 1/Pa3/s 6e-24 Glen's law fluidity coefficient (Schoof 2010)
FLUIDITY_COEFFICIENT = 6e-24
SUBGLACIAL_MELT_OPENING = 1 / ICE_DENSITY / LATENT_HEAT_FUSION
SUBGLACIAL_CREEP_PARAM = 1 * FLUIDITY_COEFFICIENT * ICE_EXPONENT ** (-ICE_EXPONENT)


def calc_ice_pressure(ice_thickness):
    return ICE_DENSITY * GRAVITY * ice_thickness

def find_nearest(array, value):
    """Finds the nearest value in an array and outputs a index.
    This function was found in 
    https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array

    Parameters:
    -----------
    array: array to be looked into
    value: single value to look for into the array

    Output:
    -------
    index of the closest value in the array
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


class MoulinShape():
    """
    This class calculates the new moulin geometry for one step


    Parameters

    head (float): initial head
    subglacial_area (float): initial subglacial cross-section area
    Qin ()

    """

    def __init__(self,                                  
                 dz = 1,
                 z_elevations = None,
                 moulin_radii = 0.2,
                 
                 temperature_profile=np.array([ZERO_KELVIN, ZERO_KELVIN]),
                 
                 ice_thickness = 500,
                 initial_head = None,
                 initial_subglacial_area = 1,
                                 
                 regional_surface_slope = 0,
                 channel_length = 15000,
                 creep_enhancement_factor = 3,

                 sigma = (-50,50),  # -50e3 #(Units??) compressive #50e3#-50e3
                 tau_xy = 0,  # -50e3#100e3 #(Units??) shear opening
                 friction_factor_OC = 0.5,
                 friction_factor_TM = 1,
                 friction_factor_SUB = 0.1,
                 fraction_pd_melting = 0.1,


                 
                 ):
        """[summary]

        Args:
            Qin ([type], optional): [description]. 
                Defaults to None.
            z_elevations (Optional, optional, list): . 
                Allows inputs of a list of floats. 
                Defaults to None.
            moulin_radii (Float, List, Tuple, optional): moulin radius in m. 
                if a float is given, moulin is assumed to be a cylinder of constant radii.
                if 2 values are given (in a list or tuple) 
                value 1 is considered the radius at the bottom 
                of the moulin and the second value is considered 
                the radius at the top of the moulin (ice surface).
                if len(moulin_radii) > 2, it must match the length of 
                z_elevations which cannot be None. 
                Defaults to 0.5 m.
            temperature_profile ([type], optional): [description]. Defaults to np.array([ZERO_KELVIN, ZERO_KELVIN]).
            tmax_in_day (int, optional): [description]. Defaults to 5.
            ice_thickness (int, optional): [description]. Defaults to 500.
            initial_head (int, optional): [description]. Defaults to 500.
            initial_subglacial_area (int, optional): [description]. Defaults to 1.
            dz (int, optional): [description]. Defaults to 1.
            dt (int, optional): [description]. Defaults to 300.
            regional_surface_slope (int, optional): [description]. Defaults to 0.
            channel_length (int, optional): [description]. Defaults to 15000.
            creep_enhancement_factor (int, optional): [description]. Defaults to 5.
            sigma (tuple, optional): [description]. Defaults to (0, 0).
            friction_factor_TM (float, optional): [description]. Defaults to 0.5.
            friction_factor_SUB (float, optional): [description]. Defaults to 0.04.
            fraction_pd_melting (float, optional): [description]. Defaults to 0.2.
            subglacial_baseflow (int, optional): [description]. Defaults to 3.
            include_ice_temperature (bool, optional): [description]. Defaults to True.
            creep (bool, optional): [description]. Defaults to True.
            elastic_deformation (bool, optional): [description]. Defaults to True.
            melt_below_head (bool, optional): [description]. Defaults to True.
            open_channel_melt (bool, optional): [description]. Defaults to False.
            potential_drop (bool, optional): [description]. Defaults to True.
            ice_motion (bool, optional): [description]. Defaults to True.
            refreezing (bool, optional): [description]. Defaults to False.
            overflow (bool, optional): [description]. Defaults to False.
        """

        
        

        #glacier ice thickness
        self.ice_thickness = ice_thickness
        #slope of the glacier surface
        self.regional_surface_slope = regional_surface_slope
        #glacier node spacing
        self.dz = dz
        #vertical array of position nodes in the moulin
        self.z = np.arange(0, self.ice_thickness + 1, self.dz)
    
        
        #Ice temperature
        self.T_ice = np.interp(self.z, np.linspace(
            0, self.z[-1], len(temperature_profile)), temperature_profile)

        # conduit length
        self.channel_length = channel_length

        # elastic deformation parameters
        self.sigma_x, self.sigma_y = sigma
        self.tau_xy = tau_xy

        # Turbulent melting parameters
        # ? may come out of __init__ later
        self.friction_factor_OC = friction_factor_OC
        self.friction_factor_TM = friction_factor_TM
        self.fraction_pd_melting = fraction_pd_melting
        self.friction_factor_SUB = friction_factor_SUB


        # creep deformation parameters
        self.creep_enhancement_factor = creep_enhancement_factor
        

        #Initializationof parameters for run1step
        #initial head and subglacial channel conduits 
        if initial_head == None:
            self.head = self.ice_thickness
        else:
            self.head = initial_head
        self.subglacial_area = initial_subglacial_area
        
        self.Qadd_total = 0
        self.dGlen = 0
        self.dGlen_cumulative = 0
        self.idx = 0
        self.listdict = []
        self.dict = defaultdict(list)


        #Initialize moulin radius
        # if no z_elevation is given, radii must either be none or less than 3 elements
        if z_elevations is None:
            if isinstance(moulin_radii, float):
                self.Mr_major = np.ones(len(self.z)) * moulin_radii
            elif len(moulin_radii) == 2:
                z_elevations = [0, self.ice_thickness]
            else:
                raise ValueError(
                    "z_elevations can not be None when inputing >2 moulin_radii values."
                )
        else:
            if len(z_elevations) != len(moulin_radii):
                raise ValueError(
                    f"Input z_elevations length {len(z_elevations)} does not match"
                    " length of moulin_radii given."
                )
        
        self.Mr_major = np.interp(
            self.z, z_elevations, moulin_radii) if z_elevations is not None else self.Mr_major
        self.Mr_minor = self.Mr_major
        self.Mx_upstream = -1 * self.Mr_major
        self.Mx_downstream = self.Mr_major
        
        
        #scalculate paramters
        self.ice_pressure = calc_ice_pressure(self.ice_thickness)
        self.ice_pressure_z = calc_ice_pressure(self.ice_thickness - self.z)
        self.C3 = (2**(5./4) / np.pi**(1/4) * np.sqrt(np.pi/(np.pi + 2))
           )/np.sqrt(WATER_DENSITY*self.friction_factor_SUB)


    def run1step(self, t, dt, meltwater_input_timeserie, 
                 subglacial_baseflow = 0,
                 head_L = None,
                 overflow=False,
                 include_ice_temperature=True,
                 creep=True,
                 elastic_deformation=True,
                 melt_below_head=True,
                 open_channel_melt=False,
                 potential_drop=True,
                 ice_motion=True,
                 refreezing=False):
        #means that you input a single value, or an array of length time
        """run the moulin shape model one timestep
    
                Args:
                    t (float or int): current time to be run. 
                    dt (float or int): timestep.
                        This is the difference between the previous timestep 
                        and the current timestep. 
                    meltwater_input_timeserie ()       
                    """    
                    

        self.head_L = head_L
        self.subglacial_baseflow = subglacial_baseflow
        self.include_ice_temperature = include_ice_temperature
        self.overflow = overflow
        self.creep = creep
        self.elastic_deformation = elastic_deformation
        self.melt_below_head = melt_below_head
        self.open_channel_melt = open_channel_melt
        self.potential_drop = potential_drop
        self.ice_motion = ice_motion
        self.refreezing = refreezing
        
        #extract current time from timeserie (s)
        self.t = t
        #extract current meltwater input from timeserie (m3/s)
        self.Qin = meltwater_input_timeserie
        #calculate lenght of current timestep (s)
        self.dt = dt#self.t[self.idx+1]-self.t[self.idx]
        


        # moulin geometry properties calculated for each time step
        ellipse_perimeter = np.pi * (3 * (self.Mr_minor + self.Mr_major) - np.sqrt(
            (3 * self.Mr_minor + self.Mr_major) * (self.Mr_minor + 3 * self.Mr_major)))
        circle_perimeter= 2 * np.pi * self.Mr_minor
        
        circle_area = np.pi * self.Mr_minor**2
        ellipse_area = np.pi * self.Mr_minor * self.Mr_major
        
        self.moulin_area = circle_area/2 + ellipse_area/2  # (m^2)        
        self.moulin_perimeter = circle_perimeter/2 + ellipse_perimeter/2  # (m)
        self.moulin_hydraulic_diameter = 4*self.moulin_area / self.moulin_perimeter
        self.moulin_hydraulic_radius = self.moulin_area/ self.moulin_perimeter
        #this is for fixes in OC
        self.circle_perimeter_OC = 2 * np.pi * self.Mr_major
        self.circle_area_OC = np.pi * self.Mr_major**2
        self.moulin_hydraulic_diameter_forOC = 4*self.circle_area_OC /self.circle_perimeter_OC
        self.moulin_hydraulic_radius_forOC = self.circle_area_OC / self.circle_perimeter_OC
        
        self.Qin_compensated = self.Qin + self.subglacial_baseflow #+ self.Qadd_total

        # Calculate head
        ################
        sol = solve_ivp(calculate_h_S_schoof,
                        # initial time and end time. We solve for one timestep.
                        [0, self.dt],
                        # initial head and channel cross-section area. Uses values in previous timestep.
                        [self.head, self.subglacial_area],
                        # args() only works with scipy>=1.4. if below, then error message: missing 5 arguments
                        args=(self.moulin_area,self.z, self.ice_thickness, self.ice_pressure,self.channel_length ,self.Qin_compensated, self.C3, self.overflow, self.head_L),
                        method='LSODA'  # solver method
                        # atol = 1e-6, #tolerance. can make it faster
                        # rtol = 1e-3,
                        # max_step = 10 #change if resolution is not enough
                        )
        # (m) moulin water level
        self.head = sol.y[0][-1]
        # (m) Channel cross-section
        self.subglacial_area = sol.y[1][-1]

        # update variables
        ##################

        # calculate meltwater discharge out of the subglacial channel
        self.Qout= self.C3*self.subglacial_area**(5/4)*np.sqrt(
            abs(WATER_DENSITY*GRAVITY*self.head/self.channel_length))
        # locate z nodes with at and under the water level
        self.wet = self.z <= self.head
        # calculate water pressure at each depth
        self.Pw_z = WATER_DENSITY*GRAVITY*(self.head-self.z)
        self.Pw_z[np.invert(self.wet)] = 0
        # calculate pressure head
        self.head_pressure = WATER_DENSITY*GRAVITY*self.head
        # calculate the pressure melting temperature
        self.Tmw = ZERO_KELVIN+0.01 - 9.8e-8 * \
            (self.head_pressure - 611.73)
        # calculate the water hydrostatic stress (OUTWARD: Positive)
        self.sigma_z = self.Pw_z - self.ice_pressure_z
        # calculate the relative friction factor. currently not active
        #friction_factor = fmm.calculate_relative_friction_factor(Mdh,Mrh,relative_roughness,type='unknown')
        # calculate head loss for melt functions
        self.dL_upstream = self.calculate_dL()
        #arrhenius ice flow low param
        self.arrhenius_ice_flow_param = self.calculate_arrhenius_param()

        # calculate moulin change
        #########################
        zeros = np.zeros(len(self.z))

        # Creep Deformation
        self.dC_major = self.calculate_creep_moulin(self.Mr_major) if self.creep else zeros
        self.dC_minor = self.calculate_creep_moulin(self.Mr_minor) if self.creep else zeros

        # Elastic deformation
        self.dE_major = self.calculate_elastic_deformation(self.Mr_major) if self.elastic_deformation else zeros
        self.dE_minor = self.calculate_elastic_deformation(self.Mr_minor) if self.elastic_deformation else zeros

        # Turbulent melting below water level  
        self.dTM = self.calculate_melt_below_head() if self.melt_below_head else zeros
        # Open channel melting
        self.dOC = self.calculate_melt_above_head_OC() if self.open_channel_melt else zeros
        # Potential drop
        self.dPD = self.calculate_melt_above_head_PD() if self.potential_drop else zeros
        # Asymmetric deformation due to Glen's Flow Law
        self.dGlen = self.calculate_iceflow_moulin() if self.ice_motion else zeros
        # Refreezing
        self.dFR = self.calculate_refreezing(self.t) if self.refreezing else zeros

        # calculate volume change
        ##########################
        self.Qadd_C = self.calculate_Q_stress_wall(self.dC_major, self.dC_minor)
        self.Qadd_E = self.calculate_Q_stress_wall(self.dE_major, self.dE_minor)
        self.Qadd_TM = self.calculate_Q_melted_wall(self.dTM)
        self.Qadd_OC = self.calculate_Q_melted_wall(self.dOC)/2
        self.Qadd_PD = self.calculate_Q_melted_wall(self.dPD)/2
        self.Qadd_total = self.Qadd_E + self.Qadd_C + self.Qadd_TM + self.Qadd_OC + self.Qadd_PD



        # calculate total radius change
        ##################################

        self.dr_major = self.dTM + self.dC_major + self.dE_major + self.dOC + self.dPD
        self.dr_minor = self.dTM + self.dC_minor + self.dE_minor + self.dPD #+ self.dOC

        # new moulin radius (USED FOR NEXT TIMESTEP)
        ###############################################

        self.Mr_major = self.Mr_major + self.dr_major
        self.Mr_minor = self.Mr_minor + self.dr_minor

        # new moulin position
        ######################
        self.Mx_upstream = self.Mx_upstream + self.dGlen - self.dr_major
        self.Mx_downstream = self.Mx_downstream + self.dGlen + self.dr_minor
        

        

        # save values in dictionnary
        self.listdict.append({'index': self.idx,
                        'delta_creep_major': self.dC_major,
                        'delta_creep_minor': self.dC_minor,
                        'delta_elastic_major': self.dE_major,
                        'delta_elastic_minor': self.dE_minor,
                        'delta_melt_below_head': self.dTM,
                        'delta_melt_above_head_open_channel': self.dOC,
                        'delta_melt_above_head_potential_drop': self.dPD,
                        'delta_ice_motion': self.dGlen,
                        'delta_refreezing':self.dFR,
                        'change_in_radius_major':self.dr_major,
                        'change_in_radius_minor':self.dr_minor,
                        'moulin_radius_major':self.Mr_major,
                        'moulin_radius_minor':self.Mr_minor,
                        'moulin_wall_position_upstream': self.Mx_upstream,
                        'moulin_wall_position_downstream': self.Mx_downstream
                        })

        self.dict['meltwater_input_moulin'].append(self.Qin)     
        self.dict['melwater_input_compensated_moulin'].append(self.Qin_compensated)
        self.dict['meltwater_output_subglacial'].append(self.Qout)
        self.dict['subglacial_cross_section_area'].append(self.subglacial_area)
        self.dict['subglacial_radius'].append(np.sqrt(self.subglacial_area*2/np.pi))
        self.dict['head'].append(self.head)
        self.dict['subglacial_baseflow'].append(self.subglacial_baseflow)
        self.dict['head_L'].append(self.head_L)
        self.dict['all_idx'].append(self.idx)
        self.dict['time'].append(self.t)
        self.time_day = np.array(self.dict['time'])/SECINDAY

        #return self.sim
        #update index for nex timestep
        self.idx = self.idx+1
        
    def plot_radius(self,axis,idx_max=-2, z_position='heq', bottom_axis=True):
        """ Plot the moulin radii timeserie for a certain position in the moulin. 
        
        Parameters:
        -----------
            axis (string) : ax .. what is that, how to describe it??
            z_position (int, string) : altitude in the moulin where the moulin radius will be display.
                                        'heq': at equilibrium head
            bottom_axis (optional, bolean) : True --> display the time axis
                                             False --> don't display time axis
        
        Example:
        --------
            fig, ax = plt.subplots()
            sim.plot_head(ax)
            """

        Mr_z = defaultdict(list)

        if z_position == 'heq':
            for idx in self.dict['all_idx']:
                idx_position = find_nearest(self.z,self.dict['head'][idx])  
                Mr_z['major'].append(self.listdict[idx]['moulin_radius_major'][idx_position])
                Mr_z['minor'].append(self.listdict[idx]['moulin_radius_minor'][idx_position])         
        else:
            idx_position = find_nearest(self.z,z_position)
            for idx in self.dict['all_idx']: 
                Mr_z['major'].append(self.listdict[idx]['moulin_radius_major'][idx_position])
                Mr_z['minor'].append(self.listdict[idx]['moulin_radius_minor'][idx_position])
            
        axis.plot(self.time_day[0,idx_max],Mr_z['major'][0,idx_max],'-',color='black') 
        axis.plot(self.time_day[0,idx_max],Mr_z['minor'][0,idx_max],'-',color='grey') 
        
        min_time = int(min(self.time_day))
        max_time = int(max(self.time_day))
        
        axis.set_ylabel('(m)',color='black')
        axis.set_xlim([min_time,max_time])
        axis.set_ylim([min(Mr_z['major']),max(Mr_z['major'])]) 
        axis.yaxis.tick_left()
        axis.yaxis.set_label_position("left")
        axis.tick_params(axis='y', labelcolor='black')
        axis.spines['top'].set_visible(False)
        
        axis.spines['right'].set_visible(False)
        axis.spines['left'].set_color('black')
        if bottom_axis == False:
            axis.spines['bottom'].set_visible(False)
            axis.axes.xaxis.set_visible(False)
        #axis.spines['left'].set_bounds()
        #axis.set_yticks(np.arange()) 
        #axis.set_yticklabels(np.arange())

        
    def plot_moulin(self, axis, idx, 
                    left_lim = -5,
                    left_bound = -4,
                    right_lim = 5,
                    right_bound = 4,
                    ground_depth=-50,
                    axis_side = 'left'
                    ):
        """ Plot moulin shape and head in the moulin for a certain timestep 
        Parameters:
        -----------
            axis
            idx (optional, int): index run to display
        """
        
        Mx_upstream = self.listdict[idx]['moulin_wall_position_upstream']
        Mx_downstream = self.listdict[idx]['moulin_wall_position_downstream']      
        
        axis.plot(Mx_upstream,self.z,color='black') #plot major axis on the left
        axis.plot(Mx_downstream,self.z,color='black')  #plot minor axis on the right

        axis.axhspan(0, self.dict['head'][idx], facecolor ='lightblue', alpha = 1,zorder=1)
        axis.axhspan(ground_depth, 0, facecolor ='peru', alpha = 1,zorder=1)
        axis.fill_betweenx(self.z,left_lim,Mx_upstream,color='aliceblue',zorder=2)
        axis.fill_betweenx(self.z,Mx_downstream,right_lim,color='aliceblue',zorder=2)      

        axis.set_ylabel('z(m)')    
        axis.set_xlabel('(m)')
        axis.set_xlim([left_lim,right_lim]) 
        axis.set_ylim([ground_depth,self.ice_thickness]) 
        axis.spines['top'].set_visible(False)
        axis.spines['bottom'].set_position(('zero')) 
        axis.spines['bottom'].set_bounds(right_bound,left_bound)
        axis.set_xticks(np.arange(left_bound,right_bound+1,2))
        axis.set_xticklabels(np.arange(left_bound,right_bound+1,2)) 
        axis.set_yticks(np.arange(0,self.ice_thickness+1,100)) 
        axis.set_yticklabels(np.arange(0,self.ice_thickness+1,100))  
        
        if axis_side == 'right':
            axis.yaxis.tick_right()
            axis.yaxis.set_label_position('right')        
            axis.spines['left'].set_visible(False)
            axis.spines['right'].set_bounds(0,self.ice_thickness)  
                      
        if axis_side == 'left':
            axis.yaxis.tick_left()
            axis.yaxis.set_label_position('left')        
            axis.spines['right'].set_visible(False)
            axis.spines['left'].set_bounds(0,self.ice_thickness)   

    def plot_head(self,axis,
                  idx_min=0,
                  idx_max=-1,
                  spine_head_min=0,
                  bottom_axis=True,
                  color='black',
                  axis_side = 'left',
                  ground_depth=-50):
        """ Plot head timeserie 
        
        Parameters:
        -----------
            axis (string) : ax .. what is that, how to describe it??
            spine_head_min (optional, int) : min value to display 
                    on the y axis on the graph
                    default value is zero (display all the way to the bed)
            bottom_axis (optional, bolean) : True --> display the time axis
                                             False --> don't display time axis
        
        Example:
        --------
            fig, ax = plt.subplots()
            MoulinShape.plot_head(ax)
            """       

        axis.plot(self.time_day[idx_min:idx_max],self.dict['head'][idx_min:idx_max],'-',color=color)        
        axis.set_ylabel('$m$',color=color)
        axis.set_xlim([self.time_day[idx_min],self.time_day[idx_max]])
        axis.set_ylim([ground_depth,self.ice_thickness]) 
        axis.yaxis.tick_left()
        axis.yaxis.set_label_position("left")
        axis.tick_params(axis='y', labelcolor=color)
        axis.spines['top'].set_visible(False)
        
        axis.set_yticks(np.arange(spine_head_min,self.ice_thickness+1,100)) 
        axis.set_yticklabels(np.arange(spine_head_min,self.ice_thickness+1,100))
        
        if axis_side == 'right':
            axis.yaxis.tick_right()
            axis.yaxis.set_label_position('right')        
            axis.spines['left'].set_visible(False)
            axis.spines['right'].set_color(color)
            axis.spines['right'].set_bounds(spine_head_min,self.ice_thickness)
            
        if axis_side == 'left':
            axis.yaxis.tick_left()
            axis.yaxis.set_label_position('left')        
            axis.spines['right'].set_visible(False)
            axis.spines['left'].set_color(color)
            axis.spines['left'].set_bounds(spine_head_min,self.ice_thickness)
        
        if bottom_axis == False:
            axis.spines['bottom'].set_visible(False)
            axis.axes.xaxis.set_visible(False)

    def plot_subglacial_radius(self,axis,
                               idx_min=0,
                               idx_max=-1,
                               bottom_axis=True,
                               color='black',
                               axis_side = 'left'):
        '''Subglacial channel'''
        axis.plot(self.time_day[idx_min:idx_max],self.dict['subglacial_radius'][idx_min:idx_max],'-',color=color)  
        axis.set_ylabel('$m$',color=color)
        axis.set_xlim([self.time_day[idx_min],self.time_day[idx_max]])
        axis.set_ylim([min(self.dict['subglacial_radius']),max(self.dict['subglacial_radius'])])
        axis.tick_params(axis='y', labelcolor=color)
        axis.spines['top'].set_visible(False)
        
        if axis_side == 'right':
            axis.yaxis.tick_right()
            axis.yaxis.set_label_position('right')        
            axis.spines['left'].set_visible(False)
            axis.spines['right'].set_color(color)
            
        if axis_side == 'left':
            axis.yaxis.tick_left()
            axis.yaxis.set_label_position('left')        
            axis.spines['right'].set_visible(False)
            axis.spines['left'].set_color(color)
            
        if bottom_axis == False:
            axis.spines['bottom'].set_visible(False)
            axis.axes.xaxis.set_visible(False)
        
    def plot_subglacial_channel(self,axis,
                                idx_min=0,
                                idx_max=-1,
                                bottom_axis=True,
                                color='black',
                                axis_side = 'left'   
                                ):
        '''Subglacial channel'''
        axis.plot(self.time_day[idx_min:idx_max],self.dict['subglacial_cross_section_area'][idx_min:idx_max],'-',color=color)  
        axis.set_ylabel('$m^2$',color=color)
        #axis.set_xlim([self.time_day[idx_min],self.time_day[idx_max]])
        axis.set_ylim([min(self.dict['subglacial_cross_section_area']),max(self.dict['subglacial_cross_section_area'])])
        axis.tick_params(axis='y', labelcolor=color)
        axis.spines['top'].set_visible(False)
        
        if axis_side == 'right':
            axis.yaxis.tick_right()
            axis.yaxis.set_label_position('right')        
            axis.spines['left'].set_visible(False)
            axis.spines['right'].set_color(color)
            
        if axis_side == 'left':
            axis.yaxis.tick_left()
            axis.yaxis.set_label_position('left')        
            axis.spines['right'].set_visible(False)
            axis.spines['left'].set_color(color)
            
        if bottom_axis == False:
            axis.spines['bottom'].set_visible(False)
            axis.axes.xaxis.set_visible(False)
        
    def plot_Qin(self,axis,
                 idx_min=0,
                 idx_max=-1,
                 bottom_axis=True,
                 color='black',
                 axis_side = 'left'
                 ):        
        '''Qin'''
        axis.plot(self.time_day[idx_min:idx_max],self.dict['meltwater_input_moulin'][idx_min:idx_max],'-',color=color)#,label='Qin') 
        axis.set_ylabel('$m^3/s$',color=color)
        axis.set_xlim([self.time_day[idx_min],self.time_day[idx_max]])
        axis.set_ylim([min(self.dict['meltwater_input_moulin']),max(self.dict['meltwater_input_moulin'])])
        axis.spines['top'].set_visible(False)
        axis.tick_params(axis='y', labelcolor=color)       
        
        if axis_side == 'right':
            axis.yaxis.tick_right()
            axis.yaxis.set_label_position('right')        
            axis.spines['left'].set_visible(False)
            axis.spines['right'].set_color(color)
            
        if axis_side == 'left':
            axis.yaxis.tick_left()
            axis.yaxis.set_label_position('left')        
            axis.spines['right'].set_visible(False)
            axis.spines['left'].set_color(color)
        if bottom_axis == False:
            axis.spines['bottom'].set_visible(False)
            axis.axes.xaxis.set_visible(False)   
            
        #ax4.set_xticks(np.round(np.linspace(min_time,20,max_time)))
        #ax4.set_xticklabels(np.round(np.linspace(min_time,20,max_time)))

    def plot_baseflow(self,axis,
                      idx_min=0,
                      idx_max=-1,
                      bottom_axis=True,
                      color='black',
                      axis_side = 'left'):        
        '''Baseflow'''
        axis.plot(self.time_day[idx_min:idx_max],self.dict['subglacial_baseflow'][idx_min:idx_max],'-',color=color)#,label='Qin') 
        axis.set_ylabel('$m^3/s$',color=color)
        axis.set_xlim([self.time_day[idx_min],self.time_day[idx_max]])
        #axis.set_ylim([min(self.dict['subglacial_baseflow']),max(self.dict['subglacial_baseflow'])])
        axis.spines['top'].set_visible(False)
        axis.tick_params(axis='y', labelcolor=color)        
        
        if axis_side == 'right':
            axis.yaxis.tick_right()
            axis.yaxis.set_label_position('right')        
            axis.spines['left'].set_visible(False)
            axis.spines['right'].set_color(color)
            
        if axis_side == 'left':
            axis.yaxis.tick_left()
            axis.yaxis.set_label_position('left')        
            axis.spines['right'].set_visible(False)
            axis.spines['left'].set_color(color)
            
        if bottom_axis == False:
            axis.spines['bottom'].set_visible(False)
            axis.axes.xaxis.set_visible(False) 

        
    def plot_Qout(self,axis,
                  idx_min=0,
                  idx_max=-1,
                  bottom_axis=True,
                  color='black',
                  axis_side = 'left'):        
        '''Qout'''
        axis.plot(self.time_day[idx_min:idx_max],self.dict['meltwater_output_subglacial'][idx_min:idx_max],'-',color=color)#,label='Qin') 
        axis.set_ylabel('$m^3/s$',color=color)        
        axis.set_xlim([self.time_day[idx_min],self.time_day[idx_max]])
        axis.set_ylim([min(self.dict['meltwater_output_subglacial']),max(self.dict['meltwater_output_subglacial'])])
        axis.tick_params(axis='y', labelcolor=color)
        axis.spines['top'].set_visible(False)
        
        if axis_side == 'right':
            axis.yaxis.tick_right()
            axis.yaxis.set_label_position('right')        
            axis.spines['left'].set_visible(False)
            axis.spines['right'].set_color(color)
            
        if axis_side == 'left':
            axis.yaxis.tick_left()
            axis.yaxis.set_label_position('left')        
            axis.spines['right'].set_visible(False)
            axis.spines['left'].set_color(color)
            
        if bottom_axis == False:
            axis.spines['bottom'].set_visible(False)
            axis.axes.xaxis.set_visible(False)   

        
    def plot_AGU_2(self,idx,t_real,h_real):
        fig = plt.figure(figsize=(10,0))
        grid = plt.GridSpec(4,4)#, wspace=-0.7)
        ax1 = fig.add_subplot(grid[0:4, 0])
        ax2 = fig.add_subplot(grid[0:4, 1:4])#, sharey=ax1)  #hw
        ax3 = fig.add_subplot(grid[2, 1:4])#ax2.twinx() #SCs
        ax4 = fig.add_subplot(grid[3, 1:4], sharex=ax2)    #Qin-Qout
        ax5 = ax4.twinx()#fig1.add_subplot(grid[2, 1:4], sharex=ax2)    #Qin-Qout
        
        self.plot_moulin(ax1,idx) 
        self.plot_head(ax2,idx_max=idx,color='steelblue',spine_head_min=200,bottom_axis=False)   
        self.plot_subglacial_channel(ax3,idx_max=idx,color='seagreen')
        self.plot_Qin(ax4,bottom_axis=False) 
        self.plot_Qout(ax5,idx_max=idx,color='orangered',bottom_axis=False) 
        #real head
        ax2.plot(t_real/3600/24,h_real,'-',color='black')        
        ax2.legend(['measured','simulated'],loc="upper left")
        #make backgroud transparent    
        ax1.patch.set_alpha(0)
        ax2.patch.set_alpha(0)
        ax3.patch.set_alpha(0)
        ax4.patch.set_alpha(0)
        ax5.patch.set_alpha(0)
        
        
            
    def plot_AGU_3(self,fig,idx,t_real,h_real):
        grid = plt.GridSpec(5,4)#, wspace=-0.7)
        ax1 = fig.add_subplot(grid[0, 1:4])#baseflow
        ax2 = fig.add_subplot(grid[1, 1:4])#Qin
        ax3 = fig.add_subplot(grid[2:5, 0])#moulin
        ax4 = fig.add_subplot(grid[2:5, 1:4])#hw
        ax5 = fig.add_subplot(grid[4:5, 1:4])#SCs

        self.plot_baseflow(ax1,idx_max=idx,color='seagreen',bottom_axis=False,axis_side = 'right')  
        self.plot_Qin(ax2,bottom_axis=False,axis_side = 'right',color='grey')  
        self.plot_moulin(ax3,idx) 
        self.plot_head(ax4,idx_max=idx,color='steelblue',spine_head_min=200,bottom_axis=False,axis_side = 'right')    
        self.plot_subglacial_radius(ax5,idx_max=idx,color='orangered',bottom_axis=True,axis_side = 'right') 

        ax4.plot(t_real/3600/24,h_real,'-',color='black') 
        
        
        l1 = ax1.legend(['Subglacial baseflow'],loc="upper right")
        l2 = ax2.legend(['Meltwater input'],loc="upper right")
        l4 = ax4.legend(['Head simulated','Head measured'],loc="upper right")    
        l5 = ax5.legend(['Subglacial radius'],loc="upper right") 
            
        for line, text in zip(l1.get_lines(), l1.get_texts()):
            text.set_color(line.get_color())
        for line, text in zip(l2.get_lines(), l2.get_texts()):
            text.set_color(line.get_color())
        for line, text in zip(l4.get_lines(), l4.get_texts()):
            text.set_color(line.get_color())
        for line, text in zip(l5.get_lines(), l5.get_texts()):
            text.set_color(line.get_color())
        
        #make backgroud transparent    
        ax1.patch.set_alpha(0)
        ax2.patch.set_alpha(0)
        ax3.patch.set_alpha(0)
        ax4.patch.set_alpha(0)
        ax5.patch.set_alpha(0)  
        
        
    def plot_AGU_4(self,fig,t_start, t_end,t_real,h_real,
                 spine_head_min=500,
                 ground_depth=-100,
                 Q_lim = [0,4],
                 SC_lim = [0,0.5],
                 display_baseflow = True
                 ):
        #fig.patch.set_facecolor('gainsboro')
        #find index based on given time
        idx_start = self.dict['time'].index(find_nearest(self.dict['time'],t_start))
        idx_end = self.dict['time'].index(find_nearest(self.dict['time'],t_end))
        t_lim = [t_start/SECINDAY,t_end/SECINDAY]
        
        #if display_baseflow == True:
        grid = plt.GridSpec(5,3)#, wspace=-0.7)
        ax1a = fig.add_subplot(grid[0, 0:2])#baseflow
        ax1b = fig.add_subplot(grid[1, 0:2])#Qin
        ax2 = fig.add_subplot(grid[2:5, 2])#moulin
        ax3 = fig.add_subplot(grid[2:5, 0:2])#hw
        ax4 = fig.add_subplot(grid[4, 0:2])#SCs
        # else:
        #     grid = plt.GridSpec(4,2)#, wspace=-0.7)
        #     ax1b = fig.add_subplot(grid[0, 0])#Qin
        #     ax2 = fig.add_subplot(grid[1:4, 1])#moulin
        #     ax3 = fig.add_subplot(grid[1:4, 0])#hw
        #     ax4 = fig.add_subplot(grid[3, 0])#SCs            
        
        #baseflow
        #if display_baseflow == True:  
        self.plot_baseflow(ax1a,
                           idx_min=idx_start,
                           idx_max=idx_end,
                           color='seagreen',
                           bottom_axis=False,
                           axis_side = 'left') 
        ax1a.set_xlim(t_lim) 
            
        #Meltwater
        self.plot_Qin(ax1b,
                       idx_min=idx_start,
                       idx_max=idx_end,
                       bottom_axis=False,
                       axis_side = 'left',
                       color='grey') 
        ax1b.set_xlim(t_lim) 
        ax1b.set_ylim(Q_lim)
        
        #Moulin
        self.plot_moulin(ax2,
                         idx_end,
                         left_lim = -11,
                         left_bound = -10,
                         right_lim = 11,
                         right_bound = 10,
                         ground_depth=ground_depth,
                         axis_side = 'right',)
        ax2.set_xticks(np.arange(-10,10+1,5))
        ax2.set_xticklabels(np.arange(-10,10+1,5)) 
        
        #Head
        self.plot_head(ax3,
                       idx_min=idx_start,
                       idx_max=idx_end,
                       color='steelblue',
                       spine_head_min=spine_head_min,
                       bottom_axis=False,
                       axis_side = 'left',
                       ground_depth=ground_depth) 
        ax3.plot(t_real/3600/24,h_real,'-',color='black') 
        ax3.set_xlim(t_lim)
        
        self.plot_subglacial_radius(ax4,
                                   idx_min=idx_start,
                                   idx_max=idx_end,
                                   color='orangered',
                                   bottom_axis=True,
                                   axis_side = 'left')         
        ax4.set_xlim(t_lim)
        ax4.set_ylim(SC_lim)

        
        #Legend
        #if display_baseflow == True:
        l1a = ax1a.legend(['Subglacial baseflow'],loc="upper left")#, bbox_to_anchor=(1, 1.5) )
        for line, text in zip(l1a.get_lines(), l1a.get_texts()):
            text.set_color(line.get_color())  
        l1b = ax1b.legend(['Meltwater input'],loc="upper left")#, bbox_to_anchor=(1, 1.5) )
        for line, text in zip(l1b.get_lines(), l1b.get_texts()):
            text.set_color(line.get_color())  
        l3 = ax3.legend(['Head simulated','Head measured'],loc="upper left")#, bbox_to_anchor=(1, 1.05) )    
        for line, text in zip(l3.get_lines(), l3.get_texts()):
            text.set_color(line.get_color())
        l4 = ax4.legend(['Subglacial radius'],loc="upper left")#, bbox_to_anchor=(1, 1.5) )
        for line, text in zip(l4.get_lines(), l4.get_texts()):
            text.set_color(line.get_color())
        ax4.patch.set_alpha(0)
        
        
        
    def plot_MGM(self,fig,t_start, t_end,
                 spine_head_min=500,
                 ground_depth=-100,
                 Q_lim = [0,4],
                 SC_lim = [0,0.5],
                 Q_fixed = False
                 ):
        fig.patch.set_facecolor('gainsboro')
        #find index based on given time
        idx_start = self.dict['time'].index(find_nearest(self.dict['time'],t_start))
        idx_end = self.dict['time'].index(find_nearest(self.dict['time'],t_end))
        t_lim = [t_start/SECINDAY,t_end/SECINDAY]
        
        if Q_fixed == True:
            grid = plt.GridSpec(3,3)#, wspace=-0.7)
            ax2 = fig.add_subplot(grid[0:3, 2])#moulin
            ax3 = fig.add_subplot(grid[0:3, 0:2])#hw
            ax4 = fig.add_subplot(grid[2, 0:2])#SCs
        else:
            grid = plt.GridSpec(4,2)#, wspace=-0.7)
            ax1 = fig.add_subplot(grid[0, 0])#Qin
            ax2 = fig.add_subplot(grid[1:4, 1])#moulin
            ax3 = fig.add_subplot(grid[1:4, 0])#hw
            ax4 = fig.add_subplot(grid[3, 0])#SCs            
        
        #Meltwater
        if Q_fixed == False:            
            self.plot_Qin(ax1,
                           idx_min=idx_start,
                           idx_max=idx_end,
                           bottom_axis=False,
                           axis_side = 'left',
                           color='grey') 
            ax1.set_xlim(t_lim) 
            ax1.set_ylim(Q_lim)
        
        #Moulin
        self.plot_moulin(ax2,
                         idx_end,
                         left_lim = -11,
                         left_bound = -10,
                         right_lim = 11,
                         right_bound = 10,
                         ground_depth=ground_depth,
                         axis_side = 'right',)
        ax2.set_xticks(np.arange(-10,10+1,5))
        ax2.set_xticklabels(np.arange(-10,10+1,5)) 
        
        #Head
        self.plot_head(ax3,
                       idx_min=idx_start,
                       idx_max=idx_end,
                       color='steelblue',
                       spine_head_min=spine_head_min,
                       bottom_axis=False,
                       axis_side = 'left',
                       ground_depth=ground_depth) 
        ax3.set_xlim(t_lim)
        
        self.plot_subglacial_radius(ax4,
                                   idx_min=idx_start,
                                   idx_max=idx_end,
                                   color='orangered',
                                   bottom_axis=True,
                                   axis_side = 'left')         
        ax4.set_xlim(t_lim)
        ax4.set_ylim(SC_lim)

        
        #Legend
        if Q_fixed == False:
            l1 = ax1.legend(['Meltwater input'],loc="upper right", bbox_to_anchor=(1, 1.5) )
            for line, text in zip(l1.get_lines(), l1.get_texts()):
                text.set_color(line.get_color())  
                ax1.patch.set_alpha(0)
        l3 = ax3.legend(['Head simulated','Head measured'],loc="upper right", bbox_to_anchor=(1, 1.05) )    
        for line, text in zip(l3.get_lines(), l3.get_texts()):
            text.set_color(line.get_color())
        l4 = ax4.legend(['Subglacial radius'],loc="upper right", bbox_to_anchor=(1, 1.5) )
        for line, text in zip(l4.get_lines(), l4.get_texts()):
            text.set_color(line.get_color())
        
        #make backgroud transparent    
        
        ax2.patch.set_alpha(0)
        ax3.patch.set_alpha(0)
        ax4.patch.set_alpha(0)
         
        
    # def plot_AGU(self,idx,t_real,h_real,
    #           spine_head_min=200
    #           ):
    
    #     time = self.time/24/3600
    #     z = self.z
    #     ice_thickness = self.ice_thickness    
    #     head = self.dict['head']
    #     subglacial_area = self.dict['subglacial_cross_section_area']   
    #     meltwater_input_moulin = self.dict['meltwater_input_moulin']
    #     #melwater_input_compensated = moulin.dict['melwater_input_compensated_moulin']
    #     meltwater_output_subglacial = self.dict['meltwater_output_subglacial'] 
    #     Mx_upstream = self.listdict[idx]['moulin_wall_position_upstream']
    #     Mx_downstream = self.listdict[idx]['moulin_wall_position_downstream']        
          
    #     fig = plt.figure(figsize=(25,8))
    #     grid = plt.GridSpec(4,4)#, wspace=-0.7)
    #     ax1 = fig.add_subplot(grid[0:4, 0])
    #     ax2 = fig.add_subplot(grid[0:4, 1:4])#, sharey=ax1)  #hw
    #     ax3 = fig.add_subplot(grid[2, 1:4])#ax2.twinx() #SCs
    #     ax4 = fig.add_subplot(grid[3, 1:4], sharex=ax2)    #Qin-Qout
    #     ax5 = ax4.twinx()#fig1.add_subplot(grid[2, 1:4], sharex=ax2)    #Qin-Qout
        
    #     # # dlim = 10
    #     # # dbound = 10
    #     # # dtick = [-10,5,0,5,10]
        
    #     time_plot=time[0:idx+1]
    #     min_time = int(min(time))
    #     max_time = int(max(time))
        
    #     ax1.clear() #clear figure content -- especially important for ploting in the loop
    #     ax2.clear()
    #     ax3.clear()
    #     ax4.clear()
    #     ax5.clear() 
        
    #     # #MOULIN SHAPE
        
    #     ax1.plot(Mx_upstream,z,color='black') #plot major axis on the left
    #     ax1.plot(Mx_downstream,z,color='black')  #plot minor axis on the right
        
    #     ax2.plot(t_real/3600/24,h_real,'-',color='black')      
    #     ax2.plot(time_plot,head[0:idx+1],'-',color='blue')   
    #     ax2.legend(['measured','simulated'],loc="upper left")
         
    #     ax3.plot(time_plot,subglacial_area[0:idx+1],'-',color='red')       
        
    #     ax4.plot(time,meltwater_input_moulin,'--',color='blue')#,label='Qin')    
    #     ax5.plot(time_plot,meltwater_output_subglacial[0:idx+1],color='red')#,'-',color='red',label='Qout')
        
    #     ax1.axhspan(0, head[idx], facecolor ='lightblue', alpha = 1,zorder=1)
    #     ax1.axhspan(-140, 0, facecolor ='peru', alpha = 1,zorder=1)
    #     ax1.fill_betweenx(z,-10,Mx_upstream,color='aliceblue',zorder=2)
    #     ax1.fill_betweenx(z,Mx_downstream,10,color='aliceblue',zorder=2)     
        
    #     '''Moulin'''
    #     ax1.set_ylabel('z(m)')    
    #     ax1.set_xlabel('(m)')
    #     ax1.set_xlim([-5,5]) 
    #     ax1.set_ylim([-50,ice_thickness]) 
    #     ax1.spines['top'].set_visible(False)
    #     ax1.spines['right'].set_visible(False)
    #     ax1.spines['bottom'].set_position(('zero')) 
    #     ax1.spines['left'].set_bounds(0,ice_thickness)
    #     ax1.spines['bottom'].set_bounds(-5,5)
    #     ax1.set_xticks([-4,-3,-2,-1,0,1,2,3,4]) 
    #     ax1.set_xticklabels([-4,-3,-2,-1,0,1,2,3,4]) 
    #     ax1.set_yticks(np.arange(0,ice_thickness+1,100)) 
    #     ax1.set_yticklabels(np.arange(0,ice_thickness+1,100))
        
    #     '''Head'''
    #     ax2.set_ylabel('Head',color='blue')
    #     ax2.set_xlim([min_time,max_time])
    #     ax2.set_ylim([-50,ice_thickness]) 
    #     ax2.yaxis.tick_left()
    #     ax2.yaxis.set_label_position("left")
    #     ax2.tick_params(axis='y', labelcolor='blue')
    #     ax2.spines['top'].set_visible(False)
    #     ax2.spines['bottom'].set_visible(False)
    #     ax2.spines['right'].set_visible(False)
    #     ax2.spines['left'].set_color('blue')
    #     ax2.axes.xaxis.set_visible(False)
    #     ax2.spines['left'].set_bounds(spine_head_min,ice_thickness)
    #     ax2.set_yticks(np.arange(spine_head_min,ice_thickness+1,100)) 
    #     ax2.set_yticklabels(np.arange(spine_head_min,ice_thickness+1,100))
        
        
    #     '''Subglacial channel'''
    #     ax3.set_ylabel('Subglacial_area',color='red')
    #     ax3.set_xlim([min_time,max_time])
    #     ax3.set_ylim([min(subglacial_area),max(subglacial_area)])
    #     ax3.yaxis.tick_right()
    #     ax3.yaxis.set_label_position("right")
    #     ax3.tick_params(axis='y', labelcolor='red')
    #     ax3.spines['top'].set_visible(False)
    #     ax3.spines['bottom'].set_visible(False)
    #     ax3.spines['left'].set_visible(False)
    #     ax3.spines['right'].set_color('red')
    #     ax3.axes.xaxis.set_visible(False)
        
        
    #     '''Qin'''
    #     ax4.set_ylabel('Qin',color='blue')
    #     ax4.set_xlim([min_time,max_time])
    #     ax4.set_ylim([min(meltwater_input_moulin),max(meltwater_input_moulin)])
    #     ax4.spines['top'].set_visible(False)
    #     ax4.spines['right'].set_visible(False)
    #     ax4.spines['left'].set_color('blue')
    #     ax4.tick_params(axis='y', labelcolor='blue')
    #     #ax4.set_xticks(np.round(np.linspace(min_time,20,max_time)))
    #     #ax4.set_xticklabels(np.round(np.linspace(min_time,20,max_time)))
        
    #     '''Qout'''
    #     ax5.set_ylabel('Qout',color='red')
    #     ax5.set_xlim([min_time,max_time])
    #     ax5.set_ylim([min(meltwater_output_subglacial),max(meltwater_output_subglacial)])
    #     ax5.yaxis.tick_right()
    #     ax5.yaxis.set_label_position("right")
    #     ax5.tick_params(axis='y', labelcolor='red')
    #     ax5.spines['top'].set_visible(False)
    #     ax5.spines['left'].set_visible(False)
    #     ax5.spines['bottom'].set_visible(False)
    #     ax5.spines['right'].set_color('red')
        
        
        
    #     #make backgroud transparent    
    #     ax1.patch.set_alpha(0)
    #     ax2.patch.set_alpha(0)
    #     ax3.patch.set_alpha(0)
    #     ax4.patch.set_alpha(0)
    #     ax5.patch.set_alpha(0)
            
        
    

    def calculate_dL(self):  # dL
        """calculates the length of wall for a defined dz"""
        dr = np.diff(self.Mx_upstream)  # !!! see with Kristin if it is okay to call it dr
        # insert value at beginning of array to match shape
        dr = np.insert(dr, 0, dr[0])
        dr[-1] = dr[-2]  # (CT) not sure why we do this one
        # create new array with the maximum of either array. In matlab, it says that: protect against negative dL
        dr = np.maximum(dr, 0)
        return np.sqrt(dr**2 + self.dz**2)
       

    def calculate_moulin_head_loss(self, Q, friction_factor, moulin_hydraulic_diameter):  # head_loss_dz
        """calculate head loss following Munson 2005"""
        # Calculates water velocity in the moulin at each node
        water_velocity = Q/self.moulin_area
        # uw[np.invert(wet)] = 0 #if there is no water in a given cell, there is no water velocity
        if (water_velocity > 9.3).any == True:  # if any value in the array is bigger than 9.3
            print('Big velocity !!! \nassigning terminal velocity of 9.3ms-1')
            # create new array with the minimum of either array. so, if uw is bigger than 9.3, then the value is replaced by 9.3
            water_velocity = np.minimum(water_velocity, 9.3)
        return((water_velocity**2) * friction_factor * self.dL_upstream) / (2 * moulin_hydraulic_diameter * GRAVITY)

    def calculate_Q_melted_wall(self, change_in_radius):
        """calculate the volume of water produced by melting of the wall
        #Notes is called Qadd_turb in MouSh """
        dA = self.moulin_perimeter * change_in_radius
        dV = ICE_DENSITY/WATER_DENSITY * np.trapz(self.z, dA)
        return dV/self.dt

    def calculate_Q_stress_wall(self, change_in_radius_major, change_in_radius_minor):
        """Calculate the change in volume """
        # # new - current
        dA_circle = np.pi * \
            (self.Mr_minor+change_in_radius_minor)**2 - np.pi * self.Mr_minor**2
        dA_ellipse = np.pi * (self.Mr_major+change_in_radius_major) * (
            self.Mr_minor+change_in_radius_minor) - np.pi * self.Mr_major * self.Mr_minor
        dA_tot = (dA_circle/2 + dA_ellipse/2)
        dV = np.trapz(self.z[self.wet], dA_tot[self.wet])
        return dV/self.dt

    # CREEP OF MOULIN WALL
    def calculate_creep_moulin(self, Mr):
        """ Calculate creep closure of a water-filled borehole  

        Parameters
        ----------
        Mr : numpy.ndarray
            The moulin radius
            It is update at each timestep with :func:`~calculate_new_moulin_wall_position`
        dt : int or float
            The time interval.
        arrhenius_ice_flow_param : float
            !!! ask Kristin to describe that
        sigma_z : numpy.ndarray
            !!! ask Kristin to describe that
        E : float
            Enhancement factor for the ice creep. 5 in matlab code.

        Returns
        -------
        dC : numpy.ndarray
            The horizontal creep closure at each vertical node.


        Notes
        -----
        Based on boreholeclosure/HomeworkProblem_Vostok3G.m which Krisin Poinar did in 2013 for crevasse model    
        Borehole 3G at Vostok by Blinov and Dmitriev (1987) and Salamatin et al (1998) from Table 4 in Talalay 
        and Hooke, 2007 (Annals)    
        boreholeclosure/HomeworkProblem_Vostok3G.m divided A by 5 in order to match measured Antarctic BH closure rates
        """
        #!!! what is epsilon_dot, ask Kristin
        # should it be 1e-3 everywhere??
        epsilon_dot = self.creep_enhancement_factor*self.arrhenius_ice_flow_param * \
            (self.sigma_z/3)**3  # this is too big. it should be 1e-3
        return Mr*np.exp(epsilon_dot*self.dt)-Mr

    # TURBULENT MELTING OF MOULIN WALL BELOW WATER LEVEL

    def calculate_melt_below_head(self):
        """
        (comment from matlab) Keep uw to a max of 9.3 m/s, artificially for now, which is the terminal velocity. 
        It was getting really large (10^50 m/s!) for areas of the moulin with near-zero cross section.
        """
        #small_radius = (self.Mr_major + self.Mr_minor) < 0.1  
        head_loss_dz_TM = self.calculate_moulin_head_loss(self.Qin, self.friction_factor_TM, self.moulin_hydraulic_diameter)

        if self.include_ice_temperature == True:
            """
            Note:
            ----
            This is modified from Jarosch & Gundmundsson (2012); Nossokoff (2013), Gulley et al. (2014), 
            Nye (1976) & Fountain & Walder (1998) to include the surrounding ice temperature """

            dM = ((WATER_DENSITY * GRAVITY * self.Qin * (head_loss_dz_TM/self.dL_upstream)) / 
                  (self.moulin_perimeter * ICE_DENSITY * (WATER_HEAT_CAPACITY * (self.Tmw - self.T_ice) + LATENT_HEAT_FUSION)))*self.dt

        else:
            """This parameterization is closer to that which is traditionally used to calculate melting within a 
            subglacial channel where the ice and water are at the same temperature"""
            dM = ((WATER_DENSITY * GRAVITY * self.Qin * (head_loss_dz_TM/self.dL_upstream)) /
                  (self.moulin_perimeter * WATER_DENSITY * LATENT_HEAT_FUSION))*self.dt
            # dM should be smoothed. use savgol or savitzky-golay or ...
        # dM[dM<0.01]=0.01
        #dM[small_radius] = 0 
        dM[~self.wet] = 0

        return dM  # [dM_major, dM_minor]

    def calculate_melt_above_head_PD(self):
        #small_radius = (self.Mr_major + self.Mr_minor) < 0.1        
        dPD = (WATER_DENSITY/ICE_DENSITY * GRAVITY/LATENT_HEAT_FUSION *
               self.Qin * self.dt / self.circle_perimeter_OC) * self.fraction_pd_melting
        #dPD[small_radius] = 0        
        dPD[self.wet] = 0

        return dPD

    def calculate_melt_above_head_OC(self):  # only for upstream!!
        # note, the friction factor can be a constante or changing in function of the hydraulic properties Mrh and Mdh
        # Lauren's way
        #small_radius = (self.Mr_major + self.Mr_minor) < 0.1 
        head_loss_dz_OC = self.calculate_moulin_head_loss(self.Qin, self.friction_factor_OC, self.moulin_hydraulic_diameter_forOC)

        remove_neg = np.zeros(len(self.dL_upstream))
        remove_neg[self.dL_upstream >= 0] = 1

        if self.include_ice_temperature == True:
            dOC = (WATER_DENSITY * GRAVITY * self.Qin * (head_loss_dz_OC / self.dL_upstream)) \
                / (self.circle_perimeter_OC * ICE_DENSITY * (WATER_HEAT_CAPACITY * (ZERO_KELVIN-self.T_ice)+LATENT_HEAT_FUSION))
        else:
            dOC = (WATER_DENSITY * GRAVITY * self.Qin * head_loss_dz_OC / self.dL_upstream) / \
                self.circle_perimeter_OC /ICE_DENSITY / LATENT_HEAT_FUSION
        #dOC[small_radius] = 0
        dOC[self.wet] = 0
        dOC_dt = dOC*self.dt*remove_neg
        return dOC_dt

    # ELASTIC DEFORMATION
    def calculate_elastic_deformation(self, Mr):
        """Calculate the horizontal deformation of the moulin due to elastic deformation. 
        Elastic deformation of a cylindrical hole in a material, with water pressure Pw and ice pressure Pi
        Based on Aadnoy 1987: Model for Fluid-Induced and In-Situ Generated Stresses in a Borehole (in rock)

        Created July 19, 2018 by Kristin Poinar for use in the moulin model 
        Modified November 1, 2019 by Kristin to fix the many errors in the equation I derived from Aadnoy.

        This solution assumes plane strain at the base of the ice sheet 
        (0 vertical strain; Aadnoy assumes the necessary vertical stress to make that true)

        Note about sigma's inputs (from matlab's code):
            - compressive (-) or extensive (+)
            - sigx-sigy is about -7x more important than tauxy in determining the surface expression of the moulin
            - we want net stress at sfc to be compressive, so sigx-sigy (+) and/or tauxy (-)
            - Strain rate data (Figure 1 of "Challenges" paper) indicate about -30 kPa mean principal stress at moulin sites


        Parameters
        ----------
        Mr_major : array
            Moulin radius in function of z (upstream)
        Mr_minor : array
            Moulin radius in function of z (downstream)
        sigma_z : float
            Vertical residual pressure
        sigma_x : float
            Horizontal stress in x direction ? (units!!!)
        sigma_y : float
            Horizontal stress in y direction ? (units!!!)
        tau_xy : float
            Shear opening (units!!!)

        Returns
        -------
        TYPE
            The change in radius for the major axis and the minor axis.
        """
        return ((1 + ICE_POISSON_RATIO)*(self.sigma_z - 0.5*(self.sigma_x + self.sigma_y)) +
                0.25 * (self.sigma_x-self.sigma_y)*(1 - 3*ICE_POISSON_RATIO - 4*ICE_POISSON_RATIO**2) +
                0.25 * self.tau_xy * (2 - 3*ICE_POISSON_RATIO - 8*ICE_POISSON_RATIO**2)) * Mr/YOUNG_ELASTIC_MODULUS


    def calculate_arrhenius_param(self):
        """Glen's Flow Law -- Cuffey and Paterson Eqn. 3.35 
        -(CT)Verify that this is true. Found info in matlab code"""
        Tfrac = 1/(self.T_ice + 7e-8*self.ice_pressure_z) - 1 / \
            (TEMPERATURE_TRANSITION + 7e-8*self.ice_pressure_z)
        Qc = LOW_CREEP_ACTIVATION_ENERGY * np.ones(len(self.T_ice))
        Qc[self.T_ice > TEMPERATURE_TRANSITION] = EFFECTIVE_CREEP_ACTIVATION_ENERGY
        return ARRHENIUS_TRANSITION * np.exp(-Qc/IDEAL_GAZ_CONSTANT * Tfrac)


    # ICE MOTION -- DEFORMATION WITH GLEN'S FLOW LAW
    def calculate_iceflow_moulin(self):
        """"Calcluate the ice motion with glen's flow low at each node z."""
        X_input = self.z
        Y_input = self.arrhenius_ice_flow_param*(self.ice_thickness-self.z)**ICE_EXPONENT
        #!!! test idea: check that the length of the cumtrapz output is the same as the other delta
        return (abs(2 * (ICE_DENSITY * GRAVITY * self.regional_surface_slope)**ICE_EXPONENT * cumtrapz(Y_input, X_input, initial=0)))*self.dt

    def calculate_refreezing(self, t):
        F = CP*(self.T_ice-ZERO_KELVIN)*2 * \
            np.sqrt(ICE_TEMPERATURE_DIFFUSION * t/np.pi)/LATENT_HEAT_FUSION
        F_prev = CP*(self.T_ice-ZERO_KELVIN)*2 * \
            np.sqrt(ICE_TEMPERATURE_DIFFUSION *
                    (t-self.dt)/np.pi)/LATENT_HEAT_FUSION
        dF = F-F_prev
        dF[~self.wet] = 0
        return dF
    
    
#doesn't work yet   
    # def plot_timeseries(self,
    #                     axis,
    #                     data,
    #                     time=None,
    #                     color='black',
    #                     ylabel=None,
    #                     y_axis_side='left',
    #                     time_axis=True,
    #                     ylim=None,
    #                     yticks=None,
    #                     tlim=None,
    #                     ):
    #     """ Plot any timeserie 
        
    #     Parameters:
    #     -----------
    #         axis (string) : ax .. what is that, how to describe it??
    #         data (string, array of lenght time) : 
    #                 'head'
    #                 'subglacial_channel'
    #                 'Qin'
    #                 'Qin_compensated'
        
    #     Example:
    #     --------
    #         fig, ax = plt.subplots()
    #         MoulinShape.plot_head(ax)
    #         """
    #     #time axis
    #     if time == None: 
    #         t = self.time_day
    #         print('t=',len(t))
    #     #x axis
    #     if data == 'head':               
    #         y = self.dict['head'] 
    #         print('head=',len(data))
    #         ylabel = 'Head ($m$)'
    #     if data == 'subglacial_channel': 
    #         y = self.dict['subglacial_cross_section_area']
    #         ylabel='Subglacial area ($m^2$)'
    #     if data == 'Qin':                
    #         y = self.dict['meltwater_input_moulin']
    #         ylabel = 'Meltwater input ($m^3/s$)'
    #     if data == 'Qin_compensated':    
    #         y = self.dict['melwater_input_compensated_moulin']
    #         ylabel = 'Total meltwater flux ($m^3/s$)'
    #     else: 
    #         y = data
    #         ylabel = ylabel
        
    #     if ylim == None: 
    #         ylim = [min(y),max(y)]
    #     if tlim == None: 
    #         tlim = [min(t),max(t)]
    #     if yticks == None: 
    #         yticks = ylim[0],ylim[1]
    
        
    #     #plot
    #     axis.plot(t,y,'-',color='blue') 
    #     axis.set_ylabel(ylabel,color=color)
    #     axis.tick_params(axis='y', labelcolor=color)
    #     axis.set_xlim(tlim)
    #     axis.set_ylim(ylim) 
    #     axis.spines['top'].set_visible(False)  
    
    #     axis.set_yticks(np.arange(yticks[0],yticks[1])) 
    #     axis.set_yticklabels(np.arange(yticks[0],yticks[1]))
    
          
    #     if y_axis_side == 'left':
    #         axis.yaxis.tick_left()
    #         axis.yaxis.set_label_position("left")
    #         axis.spines['right'].set_visible(False)
    #         axis.spines['left'].set_color(color)
            
    #     if y_axis_side == 'right': 
    #         axis.yaxis.tick_left()
    #         axis.yaxis.set_label_position("right") 
    #         axis.spines['left'].set_visible(False)
    #         axis.spines['right'].set_color(color) 
    #         axis.spines['right'].set_bounds(yticks[0])
            
    #     if time_axis == False:
    #         axis.spines['bottom'].set_visible(False)
    #         axis.axes.xaxis.set_visible(False)
        





    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
#create time array for run
def TimeStamps(time_start, time_end, timestep):     
    return np.arange(time_start, time_end+1, timestep)

def Qin_constant(time,Qin):
    """calculate constant meltwater input
    Qin_mean - constant meltwater input"""
    return np.ones(len(time))*Qin

def Qin_sinusoidal(time,Qin_mean, dQ):
    """
    dQ - amplitude of meltwater oscillation
    Qin_mean
    """
    # Qin_mean=3, dQ=0.5, period=24*3600
    return dQ*np.sin(2*PI*time / SECINDAY) + Qin_mean

def Qin_real(time, Qin_data, Qtime_data):
    return np.interp(time, Qtime_data, Qin_data)

def calculate_h_S_schoof(t, y, moulin_area, z, ice_thickness, ice_pressure, channel_length , Qin_compensated, C3, overflow, head_L):
    """Input function for ode solver to calculate moulin head and subglacial channel cross-section area.
    Equation is from Schoof 2010, subglacial channel only, without the cavities

    Parameters
    ----------
    t : array
            time
    y : 2D array
            y(0): initial head
            y(1): initial channel cross-section area
    moulin_area: array
            moulin cross-section area at each depth
    L: float
            subglacial conduit length
    Pi : array
            ice pressure
    overflow : bool #!!!how can we make this optional???
        True: enable moulin head to go above ice thickness
        False: prevent overflow of head with fixed head at 0.999H

    Returns
    -------	
    hw : array_like
            moulin head timeserie
    S : array_like
            subglacial channel cross-section area time serie  

    Notes
    -----

    References
    ----------
    .. [1] Schoof, C. Ice-sheet acceleration driven by melt supply variability. Nature 468, 803–806 (2010).

    Examples
    --------
    """

    head = y[0]  # (m) moulin head initial value
    # (m) subglacial channel cross-section area initial value
    subglacial_area = y[1]

    # sets the head to just a little smaller than the ice thickness
    if overflow == False:
        if head > ice_thickness:
            head = 0.999*ice_thickness
    # #print(head)
    # if head <= 0:
    #     head = 1

    # find moulin_area value by interpolating moulin_area value at the level of the water
    moulin_area_at_head = np.interp(head, z, moulin_area)
        
        
    if head_L == None:
        # Moulin head ODE
        dhdt = 1/moulin_area_at_head * (Qin_compensated - C3*subglacial_area**(5/4)*np.sqrt((WATER_DENSITY*GRAVITY*head)/channel_length))  
        # Channel cross-section area ODE
        melt = SUBGLACIAL_MELT_OPENING * C3 * subglacial_area**(5/4) * ((WATER_DENSITY*GRAVITY*head)/channel_length)**(3/2)
        creep = SUBGLACIAL_CREEP_PARAM * (ice_pressure - WATER_DENSITY*GRAVITY *head)**ICE_EXPONENT * subglacial_area
        dSdt =  melt-creep 
        
    else:
        N_0 = ice_pressure - WATER_DENSITY * GRAVITY * head 
        N_L = ice_pressure - WATER_DENSITY * GRAVITY * head_L
        mean_effective_pressure = (N_0 + N_L) / 2
        hydraulic_gradient = (N_L - N_0) / channel_length #WATER_DENSITY * GRAVITY * (head - head_L) / channel_length 
        Q_out = C3 * subglacial_area**(5/4) * 1/np.sqrt(abs(hydraulic_gradient)) * hydraulic_gradient
        
        # Moulin head ODE
        dhdt = (Qin_compensated - Q_out) /moulin_area_at_head
        # Channel cross-section area ODE
        dSdt = SUBGLACIAL_MELT_OPENING * abs(Q_out * hydraulic_gradient) - SUBGLACIAL_CREEP_PARAM * mean_effective_pressure**ICE_EXPONENT

    # prevents the head from getting bigger if it's close to overlow
    if overflow == False:
        if head > 0.990*ice_thickness:
            if dhdt > 0:
                dhdt = 0
    # if head <= 0:
    #     head = 1
        
    return [dhdt, dSdt]


