'''
Moulin physical model Functions for run_moulin_model.

Code developped and writen in Matlab by Kristin Poinar and Lauren Andrews, 2019-2020.
Translation in Python by Celia Trunz, 2020.
'''

import numpy as np
from scipy.integrate import cumtrapz


'''Default constants'''
T0 = 273.15 #Kelvin
rhoi = 910 #kg/m3; Ice density
rhow = 1000 #kg/m3; Water density
ki = 2.1 #J/mKs
cp = 2115 #J/kgK
Lf = 335000 #J/kg; Latent heat of fusion
g = 9.8 #m/s2; Gravity
Y = 5e9 #Pa; Young's elastic modulus (Vaughan 1995)
A = (6e-24) #1/Pa3/s; 6e-24 Glen's law fluidity coefficient (Schoof 2010)
f = 0.15 #unitless; Darcy-Weisbach friction factor (0.1 in Matt's code, 0.0375 in Schoof 2010)
n = 3 #unitless; Glen's law exponent (Schoof 2010)
# Subglacialsc model constants
c1 = 1/rhoi/Lf # units; Melt opening parameter (Schoof 2010)
c2 = 1*A*n**(-n) # units; Closure parameter (Schoof 2010)
#C.c3 = 2^(1/4) * (pi+2)^(1/2) / (pi^(1/4) * (C.rhow*C.f)^(1/2)); % units; Flux parameter (Schoof 2010)
#C.c3 = 2^(5/4) / pi**(1/4) * sqrt( pi/ ((pi+2)*rho_w*f) ) # from Matt Covington pdf, 
# for a semi-circular conduit, modified from Schoof whose equation appears to be incorrect. 
c3 = ( 2**(5./4) /np.pi**(1/4) * np.sqrt(np.pi/(np.pi + 2)) )/np.sqrt(rhow*f) # Corrected by LCA?
#C.E = 1e12;
nu = 0.3    # []; Poissons ratio for ice
#
#Glen's Flow Law -- Cuffey and Paterson Eqn. 3.35
#A = Astar * exp(-Q / R *(1/Th - 1/Tstar))
#Th = 263 + 7e-8*P
#Tstar = T + 7e-8*P

Tstar = 263
Astar = 3.5e-25
c = 7e-8
Qcless = 6e4
Qcmore = 11.5e4   
R = 8.314                # ideal gas constant
a = 7e-8
#Turbulence parameters
mu = 0.0017916 # Pa*s at 0â , https://www.engineeringtoolbox.com/water-dynamic-kinematic-viscosity-d_596.html
kw = 0.555 # J/(mKs) at 0.01â , https://www.engineeringtoolbox.com/water-liquid-gas-thermal-conductivity-temperature-pressure-d_2012.html
cw = 4210   #J / (kg * K)   heat capacity of water for unit mass,  from Jarosch & Gundmundsson (2012)
kappa = ki/rhoi/cp #this might be the temperature diffusion in ice.... confirm with Kristin and Lauren
secinday = 24*3600


#to compare with matlab inputs, arange function needs +1 in python, but not in linspace.

#INITIALISATION OF MODEL

def initiate_moulin_wall_position(Mr_major_initial,Mr_minor_initial,z):
    '''define moulin wall in relation to space
    - r_minor is the radius that controls the circle and the small axis of the elipse
    and r_major control the large axis of the elipse'''
    Mr_major = Mr_major_initial * np.ones(len(z))
    Mr_minor = Mr_minor_initial * np.ones(len(z))
    Mx_upstream = -Mr_major_initial * np.ones(len(z))
    Mx_downstream = Mr_minor_initial * np.ones(len(z))
    return [Mx_upstream,Mx_downstream, Mr_major, Mr_minor]

def generate_grid_x(dt, xmax, chebx=False):
    '''define grid resolution, spacing and length in the horizontal plane
    Default is a uniform resolution in the x axis'''
    dx = np.sqrt(kappa *dt) #(m) Diffusion lengthscale
    if chebx == True:
    #resolution increase close to the moulin -- not working yet, Ask Kristin for details
        x = np.linspace(0,np.pi/2,round(xmax/dx/10))
        x = (xmax * (1-np.sin(x)))
    if chebx == False:
    # uniform resolution in x
        x = np.arange(0,xmax+1,dx)
    else:
        print('error, please choose chebx=True or chebx=False') #(CT)do this in a better way. 
        #!! we add +1 at xmax because python consider it not inclusive. 
        #This way, entered parameters are going to match matlab input
    nx = len(x)
    dx = np.append(np.diff(x),np.diff(x)[-1]) #calculate delta x between each x and adds one.. ask kristin why.   
    return [x, nx, dx]
    
def generate_grid_z(H, dz=1):
    #z: vector of position of the points of calculations.
    #nz: number of points.
    #dz:  vertical spacing between the vertical point of calculation in the moulin
    z=np.arange(0,H+1,dz) #(m) vertical profile fo moulin point of calculation 
    nz=len(z)
        #!! we add +1 at xmax because python consider it not inclusive. 
        #This way, entered parameters are going to match matlab input
    return [z, nz, dz]

def generate_time(dt,tmax_in_day):
    tmax_in_second = tmax_in_day*24*3600
    time_vector = np.arange(dt,tmax_in_second+1,dt)#+1 is important to match matlab
    return [time_vector, tmax_in_second]
        
#def generate_vector_time(number_of_days=20,timestep_in_seconds=300):
def set_Qin(time,type, Qin_mean=3, dQ=0.5, period=24*3600):
    '''input a time array'''
    if type == 'constant':
        return Qin_mean
    if type == 'sinusoidal_celia':
        return dQ * np.sin(2*np.pi*time/period) + Qin_mean   
       
def set_ice_temperature(x,z, type='Temperate'): #T
    '''Define ice temperature
    - Create x with generate_grid_x
    - Create z with generate_grid_z
    '''   

    #define temperature in the glacier, far from the moulin
    if type == 'Temperate':
        T_far = T0 * np.ones(len(z))
    if type == 'Cool':
        Tmin = -5
        T_far = np.linspace(0,Tmin,len(z)) + T0 
        #if type == 'Cold'
    #if type == 'Luthi'
    #if type == 'HarrS2A
    #... add all the option available in matlab
    
    #define initial temperature close to moulin until the limit we set
    ones_x = np.ones(len(x)) # x direction
    T_xz = np.outer( T_far.ravel(), ones_x.ravel() ) #Ambient ice temperature everywhere to start
    T_xz[:,0]=T0 #Melting point at the moulin wall
    return [T_far,T_xz]


#SMALL FUNCTIONS

def calculate_mean_ice_temperature(T): #T_mean
    ''' mean of each row of tempearture in the x axis 
    (in matlab, it's written mean(T,2))'''
    return np.mean(T, axis=1) 

def calculate_iceflow_law_parameter(T_mean,Pi_z): #iceflow_param_glen
    '''Glen's Flow Law -- Cuffey and Paterson Eqn. 3.35 
    -(CT)Verify that this is true. Found info in matlab code'''   
    Tfrac = 1/(T_mean + a*Pi_z) - 1/(Tstar + a*Pi_z)
    Qc = Qcless * np.ones(len(T_mean))
    #Qc[T>Tstar] = Qcmore
    return Astar * np.exp(-Qc/R * Tfrac)
    
def calculate_ice_pressure_at_depth(H,z): #Pi_z
    ''' Ice pressure in function of depth
    defined as P in matlab code'''
    return rhoi*g*(H-z)


def calculate_water_pressure_at_depth(hw,z,wet): #Pw_z
    ''' Water pressure in function of depth
    Input: hw = single float of moulin hydraulic head (water level with origin at the bottom of the ice)
            z = array of float with'''                                      
    Pw_z = rhow*g*(hw-z)
    Pw_z[np.invert(wet)] = 0 # There is no stress from the water above the water level !!! here, this also puts zero in Pw_z
    return Pw_z
  
def locate_water(hw,z): #wet
    '''Boolean index of position of z that are underwater. 
    in the Matlab code, it correspond to "wet" '''
    return z <= hw

def calculate_moulin_geometry(Mx_upstream, Mx_downstream, Mr_major, Mr_minor):
    '''This is valid only for the combined demi-ellipse and demi-circle
    rlv stands for relative'''
    Diameter = Mx_downstream - Mx_upstream
    Mcs = (np.pi*Mr_minor**2)/2 + (np.pi*Mr_minor*Mr_major)/2 #(m^2) relative moulin cross-section area 
    Mpr = np.pi * (3 *(Mr_minor + Mr_major) - np.sqrt((3* Mr_minor + Mr_major) * (Mr_minor +3 * Mr_major))) #Moulin perimeter
    Mdh = (4*(np.pi * Mr_minor * Mr_major)) / Mpr #hydrualic diameter
    Mrh = (np.pi* Mr_minor * Mr_major) / Mpr # hydraulic radius
    
    #if (np.isclose(Diameter , Mr_minor - Mr_major)).any != True: #test to make sure geometry is calculated right
    #create a test with real numbers to see if geometry is calculated right
        #print('Diameter=',Diameter)
        #print('Mr_major with diameter =', Diameter-Mr_minor)
        #print('Mr_major =', Mr_major )
        #print(Mr_major[Diameter-Mr_minor - Mr_major != 0])
        #print('error in calculate_moulin_geometry \ninconsistancy with the way the two moulin radius is calculated')
    return [Mcs, Mpr, Mdh, Mrh, Diameter]
    

    
#MAIN FUNCTIONS

def calculate_Qout(S,hw,L):
    '''Discharge out of the subglacial channel'''
    return c3*S**(5/4)*np.sqrt(rhow*g*hw/L)

def calculate_water_velocity(Qout,Msc,wet):
    '''Calculates water velocity in the moulin at each node
    (is called uw in matlab code)'''
    uw = Qout/Msc
    uw[np.invert(wet)] = 0 #if there is no water in a given cell, there is no water velocity
    if (uw>9.3).any == True: # if any value in the array is bigger than 9.3
        print('Big velocity !!! \nassigning terminal velocity of 9.3ms-1')
        uw = np.minimum(uw,9.3) #create new array with the minimum of either array. so, if uw is bigger than 9.3, then the value is replaced by 9.3
    return uw

def calculate_pressure_melting_temperature(Pw_h):
    '''(Matlab comment) calculate the pressure melting temperature of ice/water 
    https://glaciers.gi.alaska.edu/sites/default/files/mccarthy/Notes_thermodyn_Aschwanden.pdf
    ''' 
    return T0+0.01 - 9.8e-8 *(Pw_h - 611.73) #!!!Ask Lauren what are those values. are they constants

def calculate_bathurst_friction_factor(Mrh,relative_roughness):
    return 10* 1/(-1.987*np.log10(relative_roughness/(5.15*Mrh)))**2

def calculate_colebrook_white_friction_factor(Mdh,relative_roughness):
    return (1/(-2*np.log10((relative_roughness/Mdh)/3.7)))**2


#WATER LEVEL
def function_subglacial_schoof(t,y,Msc,z,Pi,L,Qin):
    '''Define dimensionalized function for subglacial channel and water level in moulin
    input:
    - t = time vector
    - matrice containing initial head and channel cross-section area
    - Msc = moulin cross-section area at each depth [array]
    - L = conduit length
    - Pi = ice pressure
    output:
    - head timeserie
    - channel cross-section area time serie  
    .. code from Celia Trunz, from Schoof 2010, subglacial channel only, without the cavities'''
    hw = y[0] #(m) moulin head initial value
    SCs = y[1] #(m) subglacial channel cross-section area initial value
    Msc_hw = np.interp(hw,z,Msc)#find Msc value by interpolating Msc value at the level of the water
    dhwdt = 1/Msc_hw* ( Qin - c3*SCs**(5/4)*np.sqrt((rhow*g*hw)/L) )# Moulin head ODE
    dSCsdt = c1 * c3 * SCs**(5/4) * ((rhow*g*hw)/L)**(3/2) - c2 * ( Pi - rhow*g*hw )**n * SCs# Channel cross-section area ODE
    return [dhwdt, dSCsdt]

#CREEP
def calculate_creep_moulin(Mr_major,Mr_minor,dt,iceflow_param_glen,Pwi_z,E):
    ''' Creep closure of a water-filled borehole   
    Based on boreholeclosure/HomeworkProblem_Vostok3G.m which Krisin Poinar did in 2013 for crevasse model    
    Borehole 3G at Vostok by Blinov and Dmitriev (1987) and Salamatin et al (1998) from Table 4 in Talalay and Hooke, 2007 (Annals)    
    .. code writen in Matlab by Kristin Poinar
    '''  
    #total stress

    epsilon_dot = E*iceflow_param_glen*(Pwi_z/3)**3  #this is too big. it should be 10-3
    #boreholeclosure/HomeworkProblem_Vostok3G.m divided A by 5 in order to match measured Antarctic BH closure rates
    #Creep closure rate
    dC_major = Mr_major*np.exp(epsilon_dot*dt)-Mr_major
    dC_minor = Mr_minor*np.exp(epsilon_dot*dt)-Mr_minor
    return [dC_major,dC_minor]

def calculate_dL(Mx, dz): #dL
    '''calculates the length of wall for a defined dz'''
    dr = np.diff(Mx)#!!! see with Kristin if it is okay to call it dr
    dr = np.insert(dr,0,dr[0]) #insert value at beginning of array to match shape
    dr[-1] = dr[-2] #(CT) not sure why we do this one 
    dr = np.maximum(dr,0) #create new array with the maximum of either array. In matlab, it says that: protect against negative dL
    return np.sqrt(dr**2 + dz**2) #

def calculate_moulin_head_loss(uw, friction_factor, dL, Mdh): #head_loss_dz
    '''calculate head loss following Munson 2005'''
    return((uw**2)* friction_factor * dL) /(2 * Mdh * g)

#TURBULENT MELTING
def calculate_turbulent_melting_moulin(Mx_upstream, Mx_downstream, friction_factor, uw, Tmw, Ti, z, dz, dt, Qout, Mpr, Mdh, include_ice_temperature):
    #!!! somthing is wrong with this one! what is Mpr and should it be devided in 2???
    '''
    (comment from matlab) Keep uw to a max of 9.3 m/s, artificially for now, which is the terminal velocity. 
    It was getting really large (10^50 m/s!) for areas of the moulin with near-zero cross section.
    '''
    dL_major = calculate_dL(Mx_upstream,dz)
    dL_minor = calculate_dL(Mx_downstream,dz)
    dL = (dL_major+dL_minor)/2
    
    # head_loss_dz_major = calculate_moulin_head_loss(uw, friction_factor, dL_major, Mdh)
    # head_loss_dz_minor = calculate_moulin_head_loss(uw, friction_factor, dL_minor, Mdh)
    head_loss_dz = calculate_moulin_head_loss(uw, friction_factor, dL, Mdh)
    
    if include_ice_temperature == True:
        '''This tis modified from Jarosch & Gundmundsson (2012); Nossokoff (2013), Gulley et al. (2014), 
        Nye (1976) & Fountain & Walder (1998) to include the surrounding ice temperature '''
        dM =( (rhow * g * Qout * (head_loss_dz/dL)) / (Mpr * rhoi * (cw * (Tmw - Ti) + Lf)) )*dt
        # dM_major =( (rhow * g * Qout * (head_loss_dz_major/dL_major)) / (Mpr * rhoi * (cw * (Tmw - Ti) + Lf)) )*dt
        # dM_minor =( (rhow * g * Qout * (head_loss_dz_minor/dL_minor)) / (Mpr * rhoi * (cw * (Tmw - Ti) + Lf)) )*dt
    else :
        '''This parameterization is closer to that which is traditionally used to calculate melting within a 
        subglacial channel where the ice and water are at the same temperature'''
        dM = ( (rhow * g * Qout * (head_loss_dz/dL)) / (Mpr * rhoi * Lf) )*dt
        # dM_major = ( (rhow * g * Qout * (head_loss_dz_major/dL_major)) / (Mpr * rhoi * Lf) )*dt
        # dM_minor = ( (rhow * g * Qout * (head_loss_dz_minor/dL_minor)) / (Mpr * rhoi * Lf) )*dt
        #dM should be smoothed. use savgol or savitzky-golay or ... 
    return dM #[dM_major, dM_minor]

def calculate_volume_melted_wall(dM, z, Mpr, dt):#dM_major, dM_minor, z, Mpr, dt):
    #!!! somthing is wrong with this one! what is Mpr and should it be devided in 2???
    #calculate the volume of water produced by melting of the wall
    Vadd_turb = rhoi/rhow * np.trapz(z, Mpr/2 * dM) / dt
    # Vadd_turb_major = rhoi/rhow * np.trapz(z, Mpr/2 * dM_major) / dt
    # Vadd_turb_minor = rhoi/rhow * np.trapz(z, Mpr/2 * dM_minor) / dt
    return Vadd_turb#Vadd_turb_major + Vadd_turb_minor

#ICE MOTION -- DEFORMATION WITH GLEN'S FLOW LAW (FIND BETTER TITLE-CT)
def calculate_iceflow_moulin(Pi_z, iceflow_param_glen, regional_surface_slope, H, z, dt): #d_ice_flow
    X_input = z
    Y_input = iceflow_param_glen*(H-z)**n
    #!!! test idea: check that the length of the cumtrapz output is the same as the other delta
    return ( abs(2* (rhoi*g*regional_surface_slope)**n * cumtrapz(Y_input,X_input,initial=0) ))*dt
    

#def calculate_refreezing()

def calculate_elastic_deformation(Mr_major, Mr_minor, Pwi_z, sigma_x, sigma_y, tau_xy):
    Elastic = (1 + nu)*(Pwi_z - 0.5*(sigma_x+sigma_y)) + 0.25 \
            * (sigma_x-sigma_y)*(1 - 3*nu - 4*nu**2) + 0.25 * tau_xy * (2 - 3*nu - 8*n**2)
    dE_major = Elastic * Mr_major/Y
    dE_minor = Elastic * Mr_minor/Y
    return [dE_major, dE_minor]

def calculate_new_moulin_wall_position(Mx_upstream, Mx_downstream, dGlen_cumulative, dC=[0,0], dTM=0, dE=[0,0], dGlen=0, dOC=0):
    Mx_upstream = Mx_upstream - dC[0] - dTM- dE[0] + dGlen - dOC
    Mx_downstream = Mx_downstream + dC[1] + dTM + dE[1] + dGlen        
    if (dGlen).all == 0:
        dGlen_cumulative = 0 #this prevents from adding glen if it's been deactivated. dGlen would still be calculated in the code and unvolontary added     
    Mr_major = dGlen_cumulative - Mx_upstream #(m) relative moulin radius
    Mr_minor = Mx_downstream - dGlen_cumulative   
    return [Mx_upstream,Mx_downstream,Mr_major,Mr_minor]

def calculate_cumulative_dGlen(dGlen, dGlen_cumulative):
    return dGlen_cumulative + dGlen # dGlen[0] this seems to be unnecessary

def calculate_melt_above_head(Mr_major, Qin, dt, Mpr, wet, method='potential_drop'):
    if method=='potential_drop':
        dP = (rhow/rhoi * g/Lf * Qin * dt/Mpr)*f
        dP[wet]=0
        return dP
    
    
        

