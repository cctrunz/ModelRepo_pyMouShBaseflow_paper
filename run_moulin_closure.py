''' 
Moulin Physical Model from LC Andrews and K Poinar
Translated by Celia Trunz

Component provenance:
    - Subglacial channel: Schoof2010
    - 
'''
#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import function_moulin_model as fmm
#import plot_codes.plot_moulin_model as pmm
import pandas as pd
#import plot_codes.comprehensive_plot

secinday=24*3600

ZERO_KELVIN = 273.15
ICE_DENSITY = 910  # kg/m3; Ice density
WATER_DENSITY = 1000  # kg/m3; Water density
GRAVITY = 9.8  # m/s2; Gravity


#mts_to_cmh = 100*60*60/dt #m per timestep to mm/h : change units


def calc_overburden_pressure(ice_thickness):
    return ICE_DENSITY * GRAVITY * ice_thickness
    


class Moulin():
    """
    what this class does
    
    
    Parame
    
    """
    
    def __init__(self, 
                 hw, 
                 SCs, 
                 Qin, 
                 temperature_profile = np.array(ZERO_KELVIN, ZERO_KELVIN), 
                 Mr_major, 
                 H = 500,
                 dz=1,
                 
                 regional_surface_slope = 0,
                 L = 15000,
                 E = 5,
                 
                 sigma = (0, 0),#-50e3 #(Units??) compressive #50e3#-50e3
                 tau_xy = 0, #-50e3#100e3 #(Units??) shear opening
                 friction_factor_OC = 0.1, 
                 friction_factor_TM = 0.5, 
                 fraction_pd_melting = 0.2,
                 baseflow = 3,
                 tmax_in_day = 50,
                 dt = 300, #(s) timestep
                 ):
        self.dz = dz
        self.z = np.arange(0, self.H + 1, self.dz)
        self.T_ice = np.interp(self.z,
                               np.linspace(0, self.z[-1], len(temperature_profile)),
                               temperature_profile)


        
        
        self.hw = hw
        self.SCs = SCs
        self.Qin = Qin
        self.Mr_major = Mr_major
        self.Mr_minor = self.Mr_major
        self.moulin_position = -1 * self.Mr_major, self.Mr_major
        self.H = H
        self.overburden_pressure = calc_overburden_pressure(self.H)
        self.Pi_z = calc_overburden_pressure(self.H - self.z)
        
        
        
        
        self.regional_surface_slope = regional_surface_slope
        
        #conduit length
        self.L = L
        
        # creep enhancement factor
        self.E = E
        
        
        # elastic deformation parameters
        self.sigma_x, self.sigma_y = sigma
        self.tau_xy = tau_xy
        
        
        # Turbulent melting parameters
        #? may come out of __init__ later
        self.friction_factor_OC = friction_factor_OC
        self.friction_factor_TM = friction_factor_TM 
        self.fraction_pd_melting = fraction_pd_melting
        self.baseflow = baseflow,
        self.tmax_in_day = tmax_in_day
        self.dt = dt
        self.dGlen = 0
        self.dGlen_cumulative = 0
        self.Vadd_C = 0
        self.Vadd_E = 0
        self.Vadd_TM = 0
        
        

        
    

        
        self.iceflow_param_glen = fmm.calculate_iceflow_law_parameter(self.T_ice,self.Pi_z) #(units?) called A in matlab's code
    

        # self.time = []
    
        # -- end of __init__
    
    def run1step(self,idx):
      
        # moulin geometry properties calculated for each time step  
        Mcs, Mpr, Mdh, Mrh = fmm.calculate_moulin_geometry(self.Mr_major,self.Mr_minor)
        
        Qin_compensated = self.Qin[idx]+ self.Vadd_E + self.Vadd_C + self.baseflow
        sol = solve_ivp(fmm.calculate_h_S_schoof,
                        [0, self.dt], #initial time and end time. We solve for one timestep.
                        [hw,SCs], #initial head and channel cross-section area. Uses values in previous timestep.
                        args=(Mcs,z,Pi_H,L,Qin_compensated,H,False), #args() only works with scipy>=1.4. if below, then error message: missing 5 arguments
                        method = 'LSODA' #solver method
                        # atol = 1e-6, #tolerance. can make it faster
                        # rtol = 1e-3,
                        #max_step = 10 #change if resolution is not enough
                        )
        #(m) moulin water level
        hw = sol.y[0][-1] 
        #(m) Channel cross-section
        SCs = sol.y[1][-1] 
        
        
        Qout = fmm.calculate_Qout(SCs,hw,L)            
        wet = fmm.locate_water(hw,z) 
        Pw_z = fmm.calculate_water_pressure_at_depth(hw,z,wet) #water pressure at each depth
        Tmw = fmm.calculate_pressure_melting_temperature(Pw_z)   
        #stress_hydro = Pw_z # Water hydrostatic stress (OUTWARD: Positive)'
        sigma_z = fmm.calculate_sigma_z(Pw_z, Pi_z)
        uw_TM = fmm.calculate_water_velocity(Qout, Mcs)
        uw_OC = fmm.calculate_water_velocity(Qin[idx], Mcs)
        #friction_factor = fmm.calculate_relative_friction_factor(Mdh,Mrh,relative_roughness,type='unknown')
        dL_upstream = fmm.calculate_dL(Mx_upstream, dz)
        head_loss_dz_TM = fmm.calculate_moulin_head_loss(uw_TM, friction_factor_TM, dL_upstream, Mdh)
        head_loss_dz_OC = fmm.calculate_moulin_head_loss(uw_OC, friction_factor_OC, dL_upstream, Mdh)
        
        
        '''Calculate moulin changes for each component'''
        #Creep Deformation
        dC_major = fmm.calculate_creep_moulin(Mr_major,dt,iceflow_param_glen,sigma_z,E)
        dC_minor = fmm.calculate_creep_moulin(Mr_major,dt,iceflow_param_glen,sigma_z,E)
        Vadd_C = fmm.calculate_Q_stress_wall(dC_major,dC_minor,Mr_major,Mr_minor,z,wet,dt)    
        #Turbulent melting
        dTM = fmm.calculate_melt_below_head(dL_upstream,head_loss_dz_TM, dt, Qin[idx], Mpr, wet,include_ice_temperature=True,T_ice=T_ice,Tmw=Tmw)
        vadd_TM = fmm.calculate_Q_melted_wall(dTM, z, Mpr, dt)    
        #Refreezing
            
        #Open channel melting
        dOC = fmm.calculate_melt_above_head_OC(Mr_major,Mpr,dL_upstream,dt,head_loss_dz_OC,Qin[idx],wet,include_ice_temperature=True,T_ice=T_ice)
        vadd_OC = fmm.calculate_Q_melted_wall(dOC, z, Mpr/2, dt) 
        
        #Potential drop
        dPD = fmm.calculate_melt_above_head_PD(Mr_major, Qin[idx], dt, Mpr, wet, fraction_pd_melting)   
        vadd_PD = fmm.calculate_Q_melted_wall(dPD, z, Mpr/2, dt) 
    
        #Elastic deformation
        dE_major = fmm.calculate_elastic_deformation(Mr_major, sigma_z, sigma_x, sigma_y, tau_xy)
        dE_minor = fmm.calculate_elastic_deformation(Mr_minor, sigma_z, sigma_x, sigma_y, tau_xy)
        Vadd_E = fmm.calculate_Q_stress_wall(dE_major,dE_minor,Mr_major,Mr_minor,z,wet,dt)
        #Asymmetric deformation due to Glen's Flow Law
        dGlen = fmm.calculate_iceflow_moulin(iceflow_param_glen, regional_surface_slope, H, z, dt)
        dGlen_cumulative = fmm.calculate_cumulative_dGlen(dGlen, dGlen_cumulative)
          
        
        '''Update moulin radii'''   
        dr_major = fmm.calculate_dradius(dE=dE_major, dC=dC_major, dTM=dTM, dPD=dPD)
        dr_minor = fmm.calculate_dradius(dE=dE_minor, dC=dC_minor, dTM=dTM, dPD=dPD)
        
        
        # new moulin radius (USED FOR NEXT TIMESTEP)
        Mr_major = fmm.update_moulin_radius( Mr_major,dr_major )
        Mr_minor = fmm.update_moulin_radius( Mr_minor,dr_minor )
        # if any(Mr_major)<=0:
        #     Mr_major.all[Mr_major<=0]=0.05
        # if any(Mr_minor)<=0:
        #     Mr_minor.all[Mr_minor<=0]=0.05
        # new moulin position (USED FOR NEXT TIME STEP)
        [Mx_upstream, Mx_downstream] = fmm.update_moulin_wall_position(Mx_upstream, Mx_downstream, dr_major,dr_minor, dGlen, dGlen_cumulative)
        
        
        # self.moulin_radii.append({'index': idx, 
        #                           'moulin_radii': moulin_radii,
        #                           'Q_in': Q_in
        #                           } 
        #                          )
        # self.meltwater_flux.append({'Q_in':Q_in, 'Q_out': Q_out})
    
        '''Save values'''
        results['Mx_upstream'][idx] = Mx_upstream
        results['Mx_downstream'][idx] = Mx_downstream
        results['Mr_major'][idx] = Mr_major
        results['Mr_minor'][idx] = Mr_minor
        results['dC_major'][idx] = dC_major
        results['dC_minor'][idx] = dC_minor
        results['dTM'][idx] = dTM
        #results['dOC'][idx] = dOC
        results['dGlen'][idx] =  dGlen
        results['dGlen_cumulative'][idx] =  dGlen_cumulative
        results['dE_major'][idx] = dE_major
        results['dE_minor'][idx] = dE_minor
        results['dr_major'][idx] = dr_major
        results['dr_minor'][idx] = dr_minor
        results['Mcs'][idx] =  Mcs
        results['Mpr'][idx] = Mpr
        results['Mdh'][idx] = Mdh
        results['Mrh'][idx] = Mrh
        results['Pi_z'][idx] = Pi_z
        results['Pw_z'][idx] = Pw_z
        results['wet'][idx] = wet
        results['Tmw'][idx] = Tmw
        results['sigma_z'][idx] = sigma_z
        results['uw_TM'][idx] = uw_TM
        results['uw_OC'][idx] = uw_OC
        results['friction_factor_TM'][idx] = friction_factor_TM
        results['friction_factor_OC'][idx] = friction_factor_OC
    
        
        results['hw'][idx] = hw
        results['SCs'][idx] = SCs
        results['Qout'][idx] = Qout
        results['Qin_compensated'][idx] = Qin_compensated
        results['Vadd_E'][idx] = Vadd_E
        results['Vadd_C'][idx] = Vadd_C
        results['Vadd_TM'][idx] = Vadd_TM

