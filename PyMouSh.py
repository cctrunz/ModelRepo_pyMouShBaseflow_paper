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
from typing import Union, Optional, Tuple
#import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.integrate import cumtrapz
#import pandas as pd
#import plot_codes.comprehensive_plot

secinday = 24*3600

ZERO_KELVIN = 273.15
ICE_DENSITY = 910  # kg/m3; Ice density
WATER_DENSITY = 1000  # kg/m3; Water density
GRAVITY = 9.8  # m/s2; Gravity
# J / (kg * K)   heat capacity of water for unit mass,  from Jarosch & Gundmundsson (2012)
WATER_HEAT_CAPACITY = 4210
LATENT_HEAT_FUSION = 335000  # J/kg
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
SUBGLACIAL_MELT_OPENING = 1 / WATER_DENSITY / LATENT_HEAT_FUSION
SUBGLACIAL_CREEP_PARAM = 1 * FLUIDITY_COEFFICIENT * \
    ICE_EXPONENT ** (-ICE_EXPONENT)


def calc_overburden_pressure(ice_thickness):
    return ICE_DENSITY * GRAVITY * ice_thickness


class MoulinShape():
    """
    This class calculates the new moulin geometry for one step


    Parameters

    head (float): initial head
    subglacial_area (float): initial subglacial cross-section area
    Qin ()

    """

    def __init__(self,
                 Qin=None,
                 z_elevations: Optional(list) = None,
                 moulin_radii: Union(float, list, tuple) = 0.5,
                 temperature_profile=np.array([ZERO_KELVIN, ZERO_KELVIN]),
                 tmax_in_day=5,
                 ice_thickness=500,
                 initial_head=500,
                 initial_subglacial_area=1,
                 dz=1,
                 dt=300,

                 regional_surface_slope=0,
                 channel_length=15000,
                 creep_enhancement_factor=5,

                 sigma=(0, 0),  # -50e3 #(Units??) compressive #50e3#-50e3
                 tau_xy=0,  # -50e3#100e3 #(Units??) shear opening
                 friction_factor_OC=0.1,
                 friction_factor_TM=0.5,
                 friction_factor_SUB=0.04,
                 fraction_pd_melting=0.2,

                 baseflow=3,

                 include_ice_temperature=True,
                 creep=True,
                 elastic_deformation=True,
                 melt_below_head=True,
                 open_channel_melt=False,
                 potential_drop=True,
                 ice_motion=True,
                 refreezing=False,
                 overflow=False
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
            baseflow (int, optional): [description]. Defaults to 3.
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
        self.Qin = Qin
        self.dz = dz
        self.ice_thickness = ice_thickness
        self.head = initial_head
        self.subglacial_area = initial_subglacial_area

        self.z = np.arange(0, self.ice_thickness + 1, self.dz)
        self.T_ice = np.interp(self.z, np.linspace(
            0, self.z[-1], len(temperature_profile)), temperature_profile)

        # old
        # initial_moulin_radius = top rad, bottom radis
        # define moulin radii profile

        # if no z_elevation is given, radii must either be none or less than 3 elements
        if z_elevations is None:
            if isinstance(moulin_radii, float):
                self.Mr_moulin = np.ones(len(self.z)) * moulin_radii
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
        self.moulin_position = -1 * self.Mr_major, self.Mr_major

        self.overburden_pressure = calc_overburden_pressure(self.ice_thickness)
        self.Pi_z = calc_overburden_pressure(self.ice_thickness - self.z)

        self.tmax_in_day = tmax_in_day
        self.tmax_in_second = tmax_in_day*24*3600
        self.time = np.arange(dt, self.tmax_in_second+1, dt)

        self.regional_surface_slope = regional_surface_slope

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
        self.baseflow = baseflow,
        self.dt = dt
        self.dGlen = 0
        self.dGlen_cumulative = 0

        self.creep = creep
        self.creep_enhancement_factor = creep_enhancement_factor
        self.elastic_deformation = elastic_deformation
        self.melt_below_head = melt_below_head
        self.open_channel_melt = open_channel_melt
        self.potential_drop = potential_drop
        self.ice_motion = ice_motion
        self.refreezing = refreezing

        self.C3 = (2**(5./4) / np.pi**(1/4) * np.sqrt(np.pi/(np.pi + 2))
                   )/np.sqrt(WATER_DENSITY*self.friction_factor_SUB)

        # self.time = []

        # -- end of __init__

    def initial_moulin_radius(self,
                              moulin_radius_bottom=0.5,
                              moulin_radius_top=0.5):
        return np.linspace(moulin_radius_bottom,
                           moulin_radius_top,
                           len(self.z))

    def Qin_constant(self, Qin):
        """calculate constant meltwater input"""
        return np.ones(len(self.time))*Qin

    def Qin_sinusoidal(self, Qin_mean, dQ):
        # Qin_mean=3, dQ=0.5, period=24*3600
        return dQ * np.sin(2*np.pi*self.time/secinday) + Qin_mean

    def Qin_double_sinusoidal(self, Qin_mean, dQ):
        return dQ * np.sin(2*np.pi*self.time/secinday) + 0.1 * np.sin(np.pi*self.time/(5*secinday)) + Qin_mean

    def Qin_real(self, Qin_array, time_array):
        return np.interp(self.time, time_array, Qin_array)

    def run1step(self, idx, t):

        # moulin geometry properties calculated for each time step
        ellipse_perimeter = np.pi * (3 * (self.Mr_minor + self.Mr_major) - np.sqrt(
            (3 * self.Mr_minor + self.Mr_major) * (self.Mr_minor + 3 * self.Mr_major)))
        circle_perimeter = 2 * np.pi * self.Mr_minor
        circle_area = np.pi * self.Mr_minor**2
        ellipse_area = np.pi * self.Mr_minor * self.Mr_major
        self.moulin_area = circle_area/2 + ellipse_area/2  # (m^2)
        self.moulin_perimeter = circle_perimeter/2 + ellipse_perimeter/2  # (m)
        self.moulin_hydraulic_diameter = (
            4*(np.pi * self.Mr_minor * self.Mr_major)) / self.moulin_perimeter
        self.moulin_hydraulic_radius = (
            np.pi * self.Mr_minor * self.Mr_major) / self.moulin_perimeter

        # Calculate head
        ################

        #!!!need to make calculate_h_S_schoof in a method?? how??
        sol = solve_ivp(self.calculate_h_S_schoof,
                        # initial time and end time. We solve for one timestep.
                        [0, self.dt],
                        # initial head and channel cross-section area. Uses values in previous timestep.
                        [self.head, self.subglacial_area],
                        # args() only works with scipy>=1.4. if below, then error message: missing 5 arguments
                        args=(self),
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
        self.Qout = self.C3*self.subglacial_area**(5/4)*np.sqrt(
            WATER_DENSITY*GRAVITY*self.head/self.channel_length)
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
        self.sigma_z = self.Pw_z - self.Pi_z
        # calculate the relative friction factor. currently not active
        #friction_factor = fmm.calculate_relative_friction_factor(Mdh,Mrh,relative_roughness,type='unknown')
        # calculate head loss for melt functions
        self.dL_upstream = self.calculate_dL(self, self.Mx_upstream)

        # calculate moulin change
        #########################

        # Creep Deformation
        if self.creep == 'ON':
            self.dC_major = self.calculate_creep_moulin(self, self.Mr_major)
            self.dC_minor = self.calculate_creep_moulin(self. self.Mr_minor)
        if self.creep == 'OFF':
            self.dC_major = np.zeros(len(self.z))
            self.dC_minor = np.zeros(len(self.z))

        # Elastic deformation
        if self.creep_enhancement_factorlastic_deformation == 'ON':
            self.dE_major = self.calculate_elastic_deformation(
                self, self.Mr_major)
            self.dE_minor = self.calculate_elastic_deformation(
                self, self.Mr_minor)
        if self.creep_enhancement_factorlastic_deformation == 'OFF':
            self.dE_major = np.zeros(len(self.z))
            self.dE_minor = np.zeros(len(self.z))

        # Turbulent melting
        if self.melt_below_head == 'ON':
            self.dTM = self.calculate_melt_below_head(self, self.Qout)
        if self.melt_below_head == 'OFF':
            self.dTM = np.zeros(len(self.z))

        # Open channel melting
        if self.open_channel_melt == 'ON':
            self.dOC = self.calculate_melt_above_head_OC(self, self.Qin)
        if self.open_channel_melt == 'OFF':
            self.dOC = np.zeros(len(self.z))

        # Potential drop
        if self.potential_drop == 'ON':
            self.dPD = self.calculate_melt_above_head_PD(self, self.Qin)
        if self.potential_drop == 'OFF':
            self.dPD = np.zeros(len(self.z))

        # Asymmetric deformation due to Glen's Flow Law
        if self.ice_motion == 'ON':
            self.dGlen = self.calculate_iceflow_moulin(self)
        if self.ice_motion == 'OFF':
            self.dPD = np.zeros(len(self.z))

        # Refreezing
        if self.refreezing == 'ON':
            self.dFR = self.calculate_refreezing(self, t)
        if self.refreezing == 'OFF':
            self.dPD = np.zeros(len(self.z))

        # calculate volume change
        ##########################

        self.Vadd_C = self.calculate_Q_stress_wall(
            self, self.dC_major, self.dC_minor)
        self.Vadd_E = self.calculate_Q_stress_wall(
            self, self.dE_major, self.dE_minor)
        self.vadd_TM = self.calculate_Q_melted_wall(self, self.dTM)
        self.vadd_OC = self.calculate_Q_melted_wall(self, self.dOC)/2
        self.vadd_PD = self.calculate_Q_melted_wall(self, self.dPD)/2

        self.Qin_compensated = self.Qin[idx] + self.Vadd_E + self.Vadd_C + \
            self.vadd_TM + self.vadd_OC + self.vadd_PD + self.baseflow

        # calculate total radius change
        ##################################

        self.dr_major = self.dTM + self.dC + self.dE + self.dOC + self.dPD
        self.dr_minor = self.dTM + self.dC + self.dE + self.dOC + self.dPD

        # new moulin radius (USED FOR NEXT TIMESTEP)
        ###############################################

        self.Mr_major = self.Mr_major + self.dr_major
        self.Mr_minor = self.Mr_minor + self.dr_minor

        # new moulin position
        ######################
        self.Mx_upstream = self.Mx_upstream + self.dGlen - self.dr_major
        self.Mx_downstream = self.Mx_downstream + self.dGlen + self.dr_minor

        #

        # self.moulin_radii.append({'index': idx,
        #                           'moulin_radii': moulin_radii,
        #                           'Q_in': Q_in
        #                           }
        #                          )
        # self.meltwater_flux.append({'Q_in':Q_in, 'Q_out': Q_out})
        return self

    def calculate_dL(self, Mx):  # dL
        """calculates the length of wall for a defined dz"""
        dr = np.diff(Mx)  # !!! see with Kristin if it is okay to call it dr
        # insert value at beginning of array to match shape
        dr = np.insert(dr, 0, dr[0])
        dr[-1] = dr[-2]  # (CT) not sure why we do this one
        # create new array with the maximum of either array. In matlab, it says that: protect against negative dL
        dr = np.maximum(dr, 0)
        return np.sqrt(dr**2 + self.dz**2)

    def calculate_moulin_head_loss(self, Q, friction_factor):  # head_loss_dz
        """calculate head loss following Munson 2005"""
        # Calculates water velocity in the moulin at each node
        water_velocity = Q/self.Msc
        # uw[np.invert(wet)] = 0 #if there is no water in a given cell, there is no water velocity
        if (water_velocity > 9.3).any == True:  # if any value in the array is bigger than 9.3
            print('Big velocity !!! \nassigning terminal velocity of 9.3ms-1')
            # create new array with the minimum of either array. so, if uw is bigger than 9.3, then the value is replaced by 9.3
            water_velocity = np.minimum(water_velocity, 9.3)
        return((water_velocity**2) * friction_factor * self.dL_upstream) / (2 * self.moulin_hydraulic_diameter * GRAVITY)

    def calculate_Q_melted_wall(self, change_in_radius):
        """calculate the volume of water produced by melting of the wall
        #Notes is called Vadd_turb in MouSh """
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
        iceflow_param_glen : float
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
        epsilon_dot = self.creep_enhancement_factor*self.iceflow_param_glen * \
            (self.sigma_z/3)**3  # this is too big. it should be 1e-3
        return Mr*np.exp(epsilon_dot*self.dt)-Mr

    # TURBULENT MELTING OF MOULIN WALL BELOW WATER LEVEL

    def calculate_melt_below_head(self, Q):
        """
        (comment from matlab) Keep uw to a max of 9.3 m/s, artificially for now, which is the terminal velocity. 
        It was getting really large (10^50 m/s!) for areas of the moulin with near-zero cross section.
        """
        head_loss_dz_TM = self.calculate_moulin_head_loss(
            self, Q, self.friction_factor_TM)

        if self.include_ice_temperature == True:
            """
            Note:
            ----
            This is modified from Jarosch & Gundmundsson (2012); Nossokoff (2013), Gulley et al. (2014), 
            Nye (1976) & Fountain & Walder (1998) to include the surrounding ice temperature """

            dM = ((WATER_DENSITY * GRAVITY * Q * (head_loss_dz_TM/self.dL_upstream)) / ...
                  (self.moulin_perimeter * ICE_DENSITY * (WATER_HEAT_CAPACITY * (self.Tmw - self.T_ice) + LATENT_HEAT_FUSION)))*self.dt

        else:
            """This parameterization is closer to that which is traditionally used to calculate melting within a 
            subglacial channel where the ice and water are at the same temperature"""
            dM = ((WATER_DENSITY * GRAVITY * Q * (head_loss_dz_TM/self.dL_upstream)) /
                  (self.moulin_perimeter * WATER_DENSITY * LATENT_HEAT_FUSION))*self.dt
            # dM should be smoothed. use savgol or savitzky-golay or ...
        # dM[dM<0.01]=0.01
        dM[~self.wet] = 0

        return dM  # [dM_major, dM_minor]

    def calculate_melt_above_head_PD(self, Q):
        dPD = (WATER_DENSITY/ICE_DENSITY * GRAVITY/LATENT_HEAT_FUSION *
               self.Q * self.dt / self.Mpr) * self.fraction_pd_melting
        dPD[self.wet] = 0
        return dPD

    def calculate_melt_above_head_OC(self, Q):  # only for upstream!!
        # note, the friction factor can be a constante or changing in function of the hydraulic properties Mrh and Mdh
        # Lauren's way

        head_loss_dz_OC = self.calculate_moulin_head_loss(
            self, self.Q, self.friction_factor_OC)

        remove_neg = np.zeros(len(self.dL_upstream))
        remove_neg[self.dL_upstream >= 0] = 1

        if self.include_ice_temperature == True:
            dOC = (WATER_DENSITY * GRAVITY * Q * (head_loss_dz_OC / self.dL_upstream)) \
                / (self.moulin_perimeter * ICE_DENSITY * (WATER_HEAT_CAPACITY * (ZERO_KELVIN-self.T_ice)+LATENT_HEAT_FUSION))
        else:
            dOC = (WATER_DENSITY * GRAVITY * Q * head_loss_dz_OC / self.dL_upstream) / \
                self.moulin_perimeter/ICE_DENSITY / LATENT_HEAT_FUSION

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

    # ICE MOTION -- DEFORMATION WITH GLEN'S FLOW LAW
    def calculate_iceflow_moulin(self):
        """"Calcluate the ice motion with glen's flow low at each node z. '
        Glen's Flow Law -- Cuffey and Paterson Eqn. 3.35 
        -(CT)Verify that this is true. Found info in matlab code"""
        Tfrac = 1/(self.T_ice + 7e-8*self.Pi_z) - 1 / \
            (TEMPERATURE_TRANSITION + 7e-8*self.Pi_z)
        Qc = LOW_CREEP_ACTIVATION_ENERGY * np.ones(len(self.T_ice))
        Qc[self.T_ice > TEMPERATURE_TRANSITION] = EFFECTIVE_CREEP_ACTIVATION_ENERGY
        iceflow_param_glen = ARRHENIUS_TRANSITION * \
            np.exp(-Qc/IDEAL_GAZ_CONSTANT * Tfrac)
        X_input = self.z
        Y_input = iceflow_param_glen*(self.ice_thickness-self.z)**ICE_EXPONENT
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

    def calculate_h_S_schoof(t, y, self, Q):
        """Input function for ode solver to calculate moulin head and subglacial channel cross-section area.
        Equation is from Schoof 2010, subglacial channel only, without the cavities

        Parameters
        ----------
        t : array
                time
        y : 2D array
                y(0): initial head
                y(1): initial channel cross-section area
        Msc: array
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
        if self.overflow == False:
            if head > self.ice_thickness:
                head = 0.999*self.ice_thickness

        # find Msc value by interpolating Msc value at the level of the water
        moulin_area_at_head = np.interp(head, self.z, self.moulin_area)
        dhdt = 1/moulin_area_at_head * (self.Qin_compensated - self.C3*subglacial_area**(
            5/4)*np.sqrt((WATER_DENSITY*GRAVITY*head)/self.channel_length))  # Moulin head ODE
        dSdt = SUBGLACIAL_MELT_OPENING * self.C3 * subglacial_area**(5/4) * ((WATER_DENSITY*GRAVITY*head)/self.channel_length)**(3/2) \
            - SUBGLACIAL_CREEP_PARAM * (self.overburden_pressure - WATER_DENSITY*GRAVITY *
                                        head)**ICE_EXPONENT * subglacial_area  # Channel cross-section area ODE

        # prevents the head from getting bigger if it's close to overlow
        if self.overflow == False:
            if head > 0.990*self.ice_thickness:
                if dhdt > 0:
                    dhdt = 0

        return [dhdt, dSdt]
