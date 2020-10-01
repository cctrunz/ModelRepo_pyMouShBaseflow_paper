"""
Python implementation of the Moulin Shape (MouSh) Model
Moulin physical model Functions for run_moulin_model.

Code developped and writen in Matlab by Kristin Poinar and Lauren Andrews, 2019-2020.
Translation in Python by Celia Trunz, 2020.
"""

import numpy as np
from scipy.integrate import cumtrapz


"""Default constants"""
T0 = 273.15 #Kelvin = 0°C
rhoi = 910 #kg/m3; Ice density
rhow = 1000 #kg/m3; Water density
g = 9.8 #m/s2; Gravity
ki = 2.1 #J/mKs
cp = 2115 #J/kgK
Lf = 335000 #J/kg; Latent heat of fusion
Y = 5e9 #Pa; Young's elastic modulus (Vaughan 1995) -->shearModulus in matlab
A = 6e-24 #1/Pa3/s; 6e-24 Glen's law fluidity coefficient (Schoof 2010)
f = 0.04 #unitless; Darcy-Weisbach friction factor (0.1 in Matt's code, 0.0375 in Schoof 2010)
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



def initiate_results_dictionnary(time,z):
    """Create a new empty dictionnary to store the variable at each timestep

    Parameters
    ----------
    time: numpy.ndarray
        The timesteps generated with `generate_time`.
    z: numpy.ndarray
        The moulin nodes generated with `generate_grid_z`.

    Returns
    -------
    results: dict
        A brand new dictionnary. 
    """
    
    results={}
    results['Mx_upstream']= np.zeros([len(time),len(z)])
    results['Mx_downstream'] = np.zeros([len(time),len(z)])
    results['Mr_major']= np.zeros([len(time),len(z)])
    results['Mr_minor'] = np.zeros([len(time),len(z)])
    results['Diameter'] = np.zeros([len(time),len(z)])
    results['dC_major'] = np.zeros([len(time),len(z)]) 
    results['dC_minor'] = np.zeros([len(time),len(z)]) 
    results['dTM'] = np.zeros([len(time),len(z)]) 
    # results['dTM_major'] = np.zeros([len(time),len(z)]) 
    # results['dTM_minor'] = np.zeros([len(time),len(z)]) 
    results['dGlen'] = np.zeros([len(time),len(z)]) 
    results['dGlen_cumulative'] = np.zeros([len(time),len(z)])
    results['dE_major'] = np.zeros([len(time),len(z)]) 
    results['dE_minor'] = np.zeros([len(time),len(z)]) 
    results['dOC'] = np.zeros([len(time),len(z)]) 
    results['dPD'] = np.zeros([len(time),len(z)]) 
    results['Mcs'] = np.zeros([len(time),len(z)]) 
    results['Mpr'] = np.zeros([len(time),len(z)]) 
    results['Mdh'] = np.zeros([len(time),len(z)]) 
    results['Mrh'] = np.zeros([len(time),len(z)]) 
    results['Pi_z'] = np.zeros([len(time),len(z)]) 
    results['Pw_z'] = np.zeros([len(time),len(z)]) 
    results['wet'] = np.zeros([len(time),len(z)]) 
    results['Pw'] = np.zeros([len(time),len(z)]) 
    results['Tmw'] = np.zeros([len(time),len(z)]) 
    results['sigma_z'] = np.zeros([len(time),len(z)]) 
    results['uw'] = np.zeros([len(time),len(z)]) 
    results['friction_factor_TM'] = np.zeros([len(time),len(z)])
    results['friction_factor_OC'] = np.zeros([len(time),len(z)])


    results['hw'] = np.zeros([len(time),1])
    results['SCs'] = np.zeros([len(time),1])
    results['Qout'] = np.zeros([len(time),1])
    results['Qin_compensated'] = np.zeros([len(time),1])  
    results['Vadd_C'] = np.zeros([len(time),1])
    results['Vadd_E'] = np.zeros([len(time),1])
    results['Vadd_PD'] = np.zeros([len(time),1])
    results['Vadd_OC'] = np.zeros([len(time),1])
    results['Vadd_TM'] = np.zeros([len(time),1])

    return results

#to compare with matlab inputs, arange function needs +1 in python, but not in linspace.

#INITIALISATION OF MODEL
#########################

def generate_grid_x(dt, xmax, chebx=False):
    """Set up horizontal grid.
    Default is a uniform resolution in the x axis.

    Parameters
    ----------
    dt : int
        The time interval.
    xmax : int
        The total length of the glacier
    chebx : bool, optional
        False = uniform resolution in x
        True = resolution increase close to the moulin -- not working yet, Ask Kristin for details
    
    Returns
    -------
    x : numpy.ndarray
        The X coordinates of the glacier nodes. From 0 to xmax. At x=0 is the margin.
    nx : int/float
        The number of nodes in the x-grid.
    dx : numpy.ndarray
        The distance between each node. The distance is constant with chebx=False and variable when chebx=True.

    Example
    --------
    >>> import numpy as np
    >>> import function_moulin_model as fmm
    >>> dt = 300 #(s)
    >>> xmax = 10000 #(m)
    >>> [x, nx, dx] = fmm.generate_grid_x(dt, xmax)
    >>> x
    [array([0.00000000e+00, 1.80923255e-02, 3.61846510e-02, ...,
        1.00009490e+04, 1.00009671e+04, 1.00009852e+04])
    >>> nx
    552776
    >>> dx
    array([0.01809233, 0.01809233, 0.01809233, ..., 0.01809233, 0.01809233,
        0.01809233])]

    Notes
    -----
    The timestep controls the density of the nodes. The larger the timestep, the smaller the number of nodes. Ask Kristin Why...

    """

    dx = np.sqrt(kappa *dt) #(m) Diffusion lengthscale !!! why is dx doubled...
    if chebx == True:
        x = np.linspace(0,np.pi/2,round(xmax/dx/10))
        x = (xmax * (1-np.sin(x))) #!!!Do we need to add the +1 here?? and above?
    if chebx == False:
        x = np.arange(0,xmax+1,dx)
        #!! we add +1 at xmax because python consider it not inclusive. 
        #This way, entered parameters are going to match matlab input
    nx = len(x)
    dx = np.append(np.diff(x),np.diff(x)[-1]) #calculate delta x between each x and adds one.. ask kristin why.   
    return [x, nx, dx]
    
def generate_grid_z(H, dz=1):
    """Set up vertical grid.
    
    Parameters
    ----------
    H : float
        The ice thickness.
    dz : float, optional
        The vertical spacing between the point of calculation.

    Returns
    -------
    z : numpy.ndarray
        The Z coordinates of the moulin nodes. 0 is the top or the bottom????
    nz : float,int
        The number of nodes in the z-grid.
    dz : float
        The vertical spacing between two nodes. Distance between each node is constant.

    Example
    --------
    >>> import numpy as np
    >>> import function_moulin_model as fmm
    >>> H = 1000 #(m)
    >>> [z, nz, dz] = fmm.generate_grid_z(H)
    >>> z
    array([   0,    1,    2, ...,  998,  999, 1000])
    >>> nz
    1001
    >>> dz
    1

    """
    z=np.arange(0,H+1,dz) #(m) vertical profile fo moulin point of calculation 
    nz=len(z)
        #!! we add +1 at H because python consider it not inclusive. 
        #This way, entered parameters are going to match matlab input
    return [z, nz, dz]

def initiate_moulin_wall(z,type='linear',**kwargs):
    """Set up initial x coordinate of moulin wall nodes and initialize moulin radius.
    
    Parameters
    ----------
    Mr_major_initial : float
    Mr_minor_initial : float
    z : numpy.ndarray
    type : string
        'constant'= radius is set with Mr_major_initial and Mr_minor_initial for all z
        'linear'= radius increase or decrease linearly from bottom to top
                  REQUIRE OPTIONAL ARGUMENTS:  - Mr_top
                                               - Mr_bottom

        'custom'= radius is provided for each z position. 
                  Make sure to input an array of z values with the bottom radius at position 0 in the array.
                  REQUIRE OPTIONAL ARGUMENTS:  - Mr_major_array 
                                               - Mr_minor_array
                                            
    

    Returns
    -------
    Mx_upstream : numpy.ndarray
        The X coordinates at each moulin node, upstream of the moulin, where the stream enters. 1D array of floats.
    Mx_downstream : numpy.ndarray
        The X coordinates at each moulin node, downstream of the moulin. No stream enters here. 1D array of floats.
    Mr_major : numpy.ndarray
        The moulin radius at each node, in the larger axis, which corresponds to `Mx_upstream`. 1D array of floats.
    Mr_minor : numpy.ndarray
        The moulin radius at each node, in the smaller axis, which corresponds to `Mx_upstream`. 1D array of floats.


    Notes
    -----
    The moulin has a shape of an egg. It is composed of one ellypse and one circle. 
    r_minor is the radius that controls the circle and the small axis of the elipse
    r_major control the large axis of the elipse

    example
    -------
    >>> import numpy as np
    >>> import function_moulin_model as fmm
    >>> H = 1000 #(m)
    >>> [z, nz, dz] = fmm.generate_grid_z(H)
    >>> z
    >>> Mr_major_initial = 0
    >>> Mr_minor_initial = 0
    >>> [Mx_upstream, Mx_downstream, Mr_major, Mr_minor]= fmm.initiate_moulin_wall_position(Mr_major_initial, Mr_minor_initial,z)
    """
    
    #calculate initial moulin radius 
    # if type=='constant':
    #     Mr_major = Mr_major_initial * np.ones(len(z))
    #     Mr_minor = Mr_minor_initial * np.ones(len(z))
        
    if type=='linear':
        Mr_top = kwargs.get('Mr_top', None)
        Mr_bottom = kwargs.get('Mr_bottom', None)        
        Mr_major = np.linspace(Mr_bottom,Mr_top,len(z))
        Mr_minor = np.linspace(Mr_bottom,Mr_top,len(z))
        
    if type=='custom':
        Mr_major = kwargs.get('Mr_major_array', None)
        Mr_minor = kwargs.get('Mr_minor_array', None)

    #calculate initial moulin wall position    
    Mx_upstream = -Mr_major
    Mx_downstream = Mr_minor
        
    return [Mx_upstream,Mx_downstream, Mr_major, Mr_minor]

def generate_time(dt,tmax_in_day):
    tmax_in_second = tmax_in_day*24*3600
    time_vector = np.arange(dt,tmax_in_second+1,dt)#+1 is important to match matlab
    return [time_vector, tmax_in_second]
        
#def generate_vector_time(number_of_days=20,timestep_in_seconds=300):
def set_Qin(time,type, **kwargs ):
    """input a time array"""
    if type == 'constant':
        Qin = kwargs.get('Qin', None)
        return Qin

    if type == 'sinusoidal_celia':
        #Qin_mean=3, dQ=0.5, period=24*3600
        Qin_mean = kwargs.get('Qin_mean', None)
        dQ = kwargs.get('dQ', None)
        period = kwargs.get('period', None)
        return dQ * np.sin(2*np.pi*time/period) + Qin_mean 
    
    if type == 'double_sinusoidal':
        Qin_mean = kwargs.get('Qin_mean', None)
        dQ = kwargs.get('dQ', None)
        period = kwargs.get('period', None)
        return dQ * np.sin(2*np.pi*time/period) + 0.1 * np.sin(np.pi*time/(5*period)) + Qin_mean    
    
    if type == 'field_data':
        Qin_array = kwargs.get('Qin_array', None)
        time_array = kwargs.get('time_array', None)
        return np.interp(time,time_array,Qin_array)
     
def set_ice_temperature(x,z,ice_temperature=[T0,T0]): #T
    """Generate ice temperature in x and z direction.

    Parameters
    ----------
    x : array
        position of glacier node generated with 'generate_grid_x'
    z : array
        position of moulin nodes generated with 'generate_grid_z'
    Temperature : array (minimum 2 values)
        initial array of temperature IN KELVIN to be interpolated at each z position in the moulin.
        For field values, create a array of evenly spaced values and make sure to set the H 
        at the same value as the data, otherwise the output temperature shape only will be kept.
        T[0] --> temperature at the base
        T[-1] --> temperature at the top
        - default value are for a temperate glacier with constant ice temperature at melting point
        
        z(m)  T(K)
        1000  252
        800   258
        600   260
        400   265
        200   270
        0     273
    
        
    Returns
    -------
    T_far: 1D numpy.array
        Temperature of the glacier at specific vertical nodes far away from the moulin.
        T_far(0) is closest to the bedrock, while T_far(-1) is closest to the top of the glacier
    
    T_xz: 2D numpy.array 
    
        
    zi   | T0  T_far4   T_far4   T_far4  ---> Ice surface
    z... | T0  T_far... T_far... T_far...
    z2   | T0  T_far2   T_far2   T_far2
    z1   | T0  T_far1   T_far1   T_far1
    z0   | T0  T_far0   T_far0   T_far0  ---> Bedrock
          –––––––––––––––––––––––––––––
           x0  x1       x...     xi
           
           |                        |
           v                        v
         Moulin                  Limit of 
                              moulin influence

    """   
    #create matching z_array array Temperature vector
    z_array = np.linspace(0,z[-1],len(ice_temperature))
    #interpolate Temperature value for all the z position in the moulin
    T_far = np.interp(z,z_array,ice_temperature) #kelvin
    #define initial temperature close to moulin until the limit we set
    ones_x = np.ones(len(x)) # x direction
    T_xz = np.outer( T_far.ravel(), ones_x.ravel() ) #Ambient ice temperature everywhere to start
    T_xz[:,0]=T0 #Melting point at the moulin wall
    return [T_far,T_xz]

#SMALL FUNCTIONS
################

def calculate_mean_ice_temperature(T): #T_mean
    """ mean of each row of tempearture in the x axis 
    (in matlab, it's written mean(T,2))"""
    return np.mean(T, axis=1) 

def calculate_iceflow_law_parameter(T_mean,Pi_z): #iceflow_param_glen
    """Glen's Flow Law -- Cuffey and Paterson Eqn. 3.35 
    -(CT)Verify that this is true. Found info in matlab code"""   
    Tfrac = 1/(T_mean + a*Pi_z) - 1/(Tstar + a*Pi_z)
    Qc = Qcless * np.ones(len(T_mean))
    #Qc[T>Tstar] = Qcmore
    return Astar * np.exp(-Qc/R * Tfrac)
    
def calculate_ice_pressure_at_depth(H,z): #Pi_z
    """ Ice pressure in function of depth.
    Parameters
    ----------
    H:
    z:

    Returns
    -------
    Pi_z:

    Notes
    -----
    Pi_z is stress.cryo in MouSh matlab (this is negative in matlab)

    """
    return rhoi*g*(H-z)

def calculate_water_pressure_at_depth(hw,z,wet): #Pw_z
    """ Water pressure in function of depth.

    Parameters:
    ----------- 
    hw : float
        The moulin hydraulic head (water level with origin at the bottom of the ice).
    z :
        Nodes

    Returns
    -------
    Pw_z:

    Notes
    -----
    -Pw_z is stress.hydro in MouSh matlab 

    """  

    Pw_z = rhow*g*(hw-z)
    Pw_z[np.invert(wet)] = 0 # There is no stress from the water above the water level !!! here, this also puts zero in Pw_z
    return Pw_z

def calculate_sigma_z(Pw_z,Pi_z):
    """ Calculate total stress on the moulin in function of depth.

    Parameters
    ----------
    Pw_z : numpy.ndarray
        Water pressure in function of depth.

    Pi_z : numpy.ndarray
        Ice pressure in function of depth.

    Returns
    -------
    sigma_z : numpy.ndarray

    Notes
    -----
    Outward stress is a positive stress, example water pressure.
    Inward stress is a negative stress, example ice pressure.
    """
    return Pw_z - Pi_z

def locate_water(hw,z): #wet
    """ Identify vertical nodes that are below water level.
    Parameters
    ----------
    hw : numpy.ndarray

    z : 

    return
    ------
    wet : bool
        True: Underwater z nodes.
        False: Above water z nodes.
    
    """
    return z <= hw

def calculate_perimeter(Mr_minor,Mr_major,object_type):  
    """Calculate the perimeter of a circle or an ellipse
    
    Parameters
    ----------
    Mr_minor : float or array
        DESCRIPTION.
    Mr_major : float or array
        DESCRIPTION.
    object_type : string
        'circle': return the perimeter of a circle with Mr_minor
        'ellipse': return the perimeter of an ellipse

    Returns
    -------
    float or array
    	The perimeter.

    """
    if object_type =='ellipse':
        return np.pi * (3 *(Mr_minor + Mr_major) - np.sqrt((3* Mr_minor + Mr_major) * (Mr_minor +3 * Mr_major)))
    if object_type =='circle':
        return 2 * np.pi * Mr_minor

def calculate_area(Mr_minor,Mr_major,object_type):
    """Calculate the cross-section area of a circle or an ellipse
    Parameters
    ----------
    Mr_minor : float or array
        DESCRIPTION.
    Mr_major : float or array
        DESCRIPTION.
    object_type : string
        'circle': return the area of a circle with Mr_minor
        'ellipse': return the area of an ellipse

    Returns
    -------
    float or array
        The area.
    """
    
    if object_type =='ellipse':
        return np.pi * Mr_minor * Mr_major
    if object_type =='circle':
        return np.pi * Mr_minor**2

def calculate_moulin_geometry(Mx_upstream, Mx_downstream, Mr_major, Mr_minor):
    """ Calculate the moulin size for each node z.

    Parameters
    ----------
    Mx_upstream : numpy.ndarray
        The x coordinate of the moulin wall in the upstream direction.
    Mx_downstream : numpy.ndarray
        The x coordinate of the moulin wall in the downstream direction.
    Mr_major : numpy.ndarray
        The larger moulin radius (upstream direction).
    Mr_minor : numpy.ndarray
        The smaller moulin radius (downstream direction).

    Returns
    -------
    Mcs : numpy.ndarray
        The moulin moulin cross-section for an ellipse
    Mpr : numpy.ndarray
        The moulin perimeter.
    Mdh : numpy.ndarray
        The moulin hydraulic diameter
    Mrh : numpy.ndarray
        The moulin hydraulic radius
    Diameter : numpy.ndarray
        The long diameter. This is a control diameter to make sure that we are not makin errors

    Notes
    -----
    Ramanujan formula to approximate the perimeter of an ellipse.
    ... say where the Mdh and Mrh come from
    This is valid only for the combined demi-ellipse and demi-circle.

    """
    
    circle_area = calculate_area(Mr_minor,Mr_major,'circle')
    circle_perimeter = calculate_perimeter(Mr_minor,Mr_major,'circle')
    ellipse_area = calculate_area(Mr_minor,Mr_major,'ellipse')
    ellipse_perimeter = calculate_perimeter(Mr_minor,Mr_major,'ellipse')
    
    Mcs = circle_area/2 + ellipse_area/2 #(m^2)
    Mpr = circle_perimeter/2 + ellipse_perimeter/2 #(m)
    Mdh = (4*(np.pi * Mr_minor * Mr_major)) / Mpr #
    Mrh = (np.pi* Mr_minor * Mr_major) / Mpr #
    Diameter = Mx_downstream - Mx_upstream
    
    #if (np.isclose(Diameter , Mr_minor - Mr_major)).any != True: #test to make sure geometry is calculated right
    #create a test with real numbers to see if geometry is calculated right
        #print('Diameter=',Diameter)
        #print('Mr_major with diameter =', Diameter-Mr_minor)
        #print('Mr_major =', Mr_major )
        #print(Mr_major[Diameter-Mr_minor - Mr_major != 0])
        #print('error in calculate_moulin_geometry \ninconsistancy with the way the two moulin radius is calculated')
    return [Mcs, Mpr, Mdh, Mrh, Diameter]



def calculate_water_velocity(Qout,Msc,wet):
    """Calculates water velocity in the moulin at each node
    (is called uw in matlab code)"""
    uw = Qout/Msc
    uw[np.invert(wet)] = 0 #if there is no water in a given cell, there is no water velocity
    if (uw>9.3).any == True: # if any value in the array is bigger than 9.3
        print('Big velocity !!! \nassigning terminal velocity of 9.3ms-1')
        uw = np.minimum(uw,9.3) #create new array with the minimum of either array. so, if uw is bigger than 9.3, then the value is replaced by 9.3
    return uw

def calculate_pressure_melting_temperature(Pw_h):
    """(Matlab comment) calculate the pressure melting temperature of ice/water 
    https://glaciers.gi.alaska.edu/sites/default/files/mccarthy/Notes_thermodyn_Aschwanden.pdf
    """ 
    return T0+0.01 - 9.8e-8 *(Pw_h - 611.73) #!!!Ask Lauren what are those values. are they constants


def calculate_relative_friction_factor(Mdh,Mrh,relative_roughness,type='unknown'): #fR
    if type=='unknown':
        return 1/(-2*np.log10(relative_roughness/(5.15*Mrh)))**2
    if type=='bathurst':
        return 10* 1/(-1.987*np.log10(relative_roughness/(5.15*Mrh)))**2
    if type=='colebrook_white':
        return (1/(-2*np.log10((relative_roughness/Mdh)/3.7)))**2

def calculate_alpha(H):
    """Calculate alpha (the regional surface slope) in function of H"""
    tau = 100e3 
    return tau/rhoi/g/H

def calculate_L(H):
    """Calculate the channel length in function of H for a idealized ice sheet profile. 
    from the function makeicesheetgeom.m in MouSh. 
    """
    # !!! Find more references.
    # !!! add other type of profiles?? example the square root one?
    tau = 100e3 #don't know what this is. ask Kristin and Lauren
    return H**2 * g * rhoi /2 /tau

def calculate_h_S_schoof(t,y,Msc,z,Pi,L,Qin,H,overflow):
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

    >>>sol = solve_ivp(fmm.calculate_h_S_schoof,
    ...                [0, dt], #initial time and end time. We solve for one timestep.
    ...                [hw,SCs], #initial head and channel cross-section area. Uses values in previous timestep.
    ...                args = (Mcs,z,Pi_H,L,Qin[idx]), #args() only works with scipy>=1.4. if older, then error message: missing 5 arguments
    ...                method = 'LSODA' #solver method
    ...                # atol = 1e-6, #tolerance. can make it faster
    ...                # rtol = 1e-3,
    ...                #max_step = 10 #change if resolution is not enough
    ...                )
    hw = sol.y[0][-1]  #(m) moulin water level
    SCs = sol.y[1][-1] #(m) Channel cross-section

    """

    hw = y[0] #(m) moulin head initial value
    SCs = y[1] #(m) subglacial channel cross-section area initial value

    #sets the head to just a little smaller than the ice thickness
    if overflow==False:
        if hw>H:
            hw=0.999*H

    Msc_hw = np.interp(hw,z,Msc)#find Msc value by interpolating Msc value at the level of the water
    dhwdt = 1/Msc_hw* ( Qin - c3*SCs**(5/4)*np.sqrt((rhow*g*hw)/L) )# Moulin head ODE
    dSCsdt = c1 * c3 * SCs**(5/4) * ((rhow*g*hw)/L)**(3/2) - c2 * ( Pi - rhow*g*hw )**n * SCs# Channel cross-section area ODE

    #prevents the head from getting bigger if it's close to overlow
    if overflow==False:
        if hw>0.990*H:
            if dhwdt>0:
                dhwdt=0

    return [dhwdt, dSCsdt]


# MOULIN SHAPE MODULES
#######################

#CREEP OF MOULIN WALL
def calculate_creep_moulin(Mr_major,Mr_minor,dt,iceflow_param_glen,sigma_z,E):
    """ Calculate creep closure of a water-filled borehole  
    
    Parameters
    ----------
    Mr_major : numpy.ndarray
        The moulin long radius (upstream) where the stream enters the moulin. 
        It is update at each timestep with :func:`~calculate_new_moulin_wall_position`
        The moulin short radius (downstream) where the stream enters the moulin
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
    dC_major : numpy.ndarray
        The horizontal creep closure at each vertical node for the upstream radius.
    dC_minor : numpy.ndarray
        The horizontal creep closure at each vertical node for the downstream radius.

    Notes
    -----
    Based on boreholeclosure/HomeworkProblem_Vostok3G.m which Krisin Poinar did in 2013 for crevasse model    
    Borehole 3G at Vostok by Blinov and Dmitriev (1987) and Salamatin et al (1998) from Table 4 in Talalay 
    and Hooke, 2007 (Annals)    
    boreholeclosure/HomeworkProblem_Vostok3G.m divided A by 5 in order to match measured Antarctic BH closure rates
    """  
    #!!! what is epsilon_dot, ask Kristin
    # should it be 1e-3 everywhere??
    epsilon_dot = E*iceflow_param_glen*(sigma_z/3)**3  #this is too big. it should be 1e-3   
    dC_major = Mr_major*np.exp(epsilon_dot*dt)-Mr_major
    dC_minor = Mr_minor*np.exp(epsilon_dot*dt)-Mr_minor
    return [dC_major,dC_minor]

def calculate_dL(Mx, dz): #dL
    """calculates the length of wall for a defined dz"""
    dr = np.diff(Mx)#!!! see with Kristin if it is okay to call it dr
    dr = np.insert(dr,0,dr[0]) #insert value at beginning of array to match shape
    dr[-1] = dr[-2] #(CT) not sure why we do this one 
    dr = np.maximum(dr,0) #create new array with the maximum of either array. In matlab, it says that: protect against negative dL
    return np.sqrt(dr**2 + dz**2) #

def calculate_moulin_head_loss(uw, friction_factor, dL, Mdh): #head_loss_dz
    """calculate head loss following Munson 2005"""
    return((uw**2)* friction_factor * dL) /(2 * Mdh * g)

#TURBULENT MELTING OF MOULIN WALL BELOW WATER LEVEL
def calculate_melt_below_head(Mx_upstream, Mx_downstream, friction_factor, uw, Tmw, z, dz, dt, Qout, Mpr, Mdh,wet,**kwargs):
    #!!! somthing is wrong with this one! what is Mpr and should it be devided in 2???
    """
    (comment from matlab) Keep uw to a max of 9.3 m/s, artificially for now, which is the terminal velocity. 
    It was getting really large (10^50 m/s!) for areas of the moulin with near-zero cross section.
    """
    dL_major = calculate_dL(Mx_upstream,dz)
    dL_minor = calculate_dL(Mx_downstream,dz)
    dL = (dL_major+dL_minor)/2
    
    # head_loss_dz_major = calculate_moulin_head_loss(uw, friction_factor, dL_major, Mdh)
    # head_loss_dz_minor = calculate_moulin_head_loss(uw, friction_factor, dL_minor, Mdh)
    head_loss_dz = calculate_moulin_head_loss(uw, friction_factor, dL, Mdh)
    
    include_ice_temperature = kwargs.get('include_ice_temperature', None)
    if include_ice_temperature == True:
        """
        Note:
        ----
        This is modified from Jarosch & Gundmundsson (2012); Nossokoff (2013), Gulley et al. (2014), 
        Nye (1976) & Fountain & Walder (1998) to include the surrounding ice temperature """
        T_far = kwargs.get('T_far', None)
        dM =( (rhow * g * Qout * (head_loss_dz/dL)) / (Mpr * rhoi * (cw * (Tmw - T_far) + Lf)) )*dt
        # dM_major =( (rhow * g * Qout * (head_loss_dz_major/dL_major)) / (Mpr * rhoi * (cw * (Tmw - Ti) + Lf)) )*dt
        # dM_minor =( (rhow * g * Qout * (head_loss_dz_minor/dL_minor)) / (Mpr * rhoi * (cw * (Tmw - Ti) + Lf)) )*dt
    else :
        """This parameterization is closer to that which is traditionally used to calculate melting within a 
        subglacial channel where the ice and water are at the same temperature"""
        dM = ( (rhow * g * Qout * (head_loss_dz/dL)) / (Mpr * rhoi * Lf) )*dt
        # dM_major = ( (rhow * g * Qout * (head_loss_dz_major/dL_major)) / (Mpr * rhoi * Lf) )*dt
        # dM_minor = ( (rhow * g * Qout * (head_loss_dz_minor/dL_minor)) / (Mpr * rhoi * Lf) )*dt
        #dM should be smoothed. use savgol or savitzky-golay or ... 
        dM[~wet]=0
    return dM #[dM_major, dM_minor]

def calculate_melt_above_head_PD(Mr_major, Qin, dt, Mpr, wet, fraction_pd_melting):
        dPD = (rhow/rhoi * g/Lf * Qin * dt / Mpr) * fraction_pd_melting
        dPD[wet]=0
        return dPD

def calculate_melt_above_head_OC(Mr_major,Mx_upstream,dz,friction_factor,Qin,wet,**kwargs): #only for upstream!!
    #note, the friction factor can be a constante or changing in function of the hydraulic properties Mrh and Mdh    
    dL_major = calculate_dL(Mx_upstream,dz)
    #Lauren's way
    area_ell = np.pi * Mr_major**2
    Mp = 2*np.pi*Mr_major
    Dh = (4*area_ell)/Mp
    Rh = area_ell/Mp
    #expected headloss based on the discharge
    hL = (Qin/area_ell)**2 * friction_factor * dL_major /2 /Dh /g
    
    remove_neg = np.zeros(len(dL_major))
    remove_neg[dL_major>=0]=1    
        
    include_ice_temperature = kwargs.get('include_ice_temperature', None)
    if include_ice_temperature==True:
        T_far = kwargs.get('T_far', None)
        dOC_dt = rhow * g * Qin * hL /dL_major /Mp /rhoi /cw /(T0-T_far+Lf)
    else:
        dOC_dt = (rhow * g * Qin * hL /dL_major) /Mp /rhoi /Lf
        
    dOC_dt[wet]=0
    dOC_dt=dOC_dt*remove_neg
    return dOC_dt
        



def calculate_Q_melted_wall(dmelt, z, Mpr, dt):
    #calculate the volume of water produced by melting of the wall
    #Notes is called Vadd_turb in MouSh 
    dA = Mpr * dmelt
    dV = rhoi/rhow * np.trapz(z, dA) 
    return dV/dt



def calculate_Q_stress_wall(dstress_major,dstress_minor,Mr_major,Mr_minor,z,wet,dt):
	dA_circle = np.pi * (Mr_minor+dstress_minor)**2 - np.pi * Mr_minor**2 # new - current
	dA_ellipse = np.pi * (Mr_major+dstress_major) * (Mr_minor+dstress_minor) - np.pi * Mr_major * Mr_minor # new - current
	dA_tot = (dA_circle/2 + dA_ellipse/2)
	dV = np.trapz(z[wet], dA_tot[wet])
	return dV/dt

#ICE MOTION -- DEFORMATION WITH GLEN'S FLOW LAW 
def calculate_iceflow_moulin(iceflow_param_glen, regional_surface_slope, H, z, dt): #d_ice_flow
    X_input = z
    Y_input = iceflow_param_glen*(H-z)**n
    #!!! test idea: check that the length of the cumtrapz output is the same as the other delta
    return ( abs(2* (rhoi*g*regional_surface_slope)**n * cumtrapz(Y_input,X_input,initial=0) ))*dt
    
#def calculate_refreezing()

#ELASTIC DEFORMATION
def calculate_elastic_deformation(Mr_major, Mr_minor, sigma_z, sigma_x, sigma_y, tau_xy):
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
    Elastic = (1 + nu)*(sigma_z - 0.5*(sigma_x+sigma_y)) + 0.25 * (sigma_x-sigma_y)*(1 - 3*nu - 4*nu**2) + 0.25 * tau_xy * (2 - 3*nu - 8*nu**2)
    dE_major = Elastic * Mr_major/Y
    dE_minor = Elastic * Mr_minor/Y
    return [dE_major, dE_minor]

def calculate_dradius(dC=[0,0], dTM=0, dE=[0,0], dOC=0, dPD=0):
    dr_major = dTM + dC[0] + dE[0] + dOC #upstream
    dr_minor = dTM + dC[1] + dE[1] + dPD #downstream
    return [dr_major,dr_minor]

#MOULIN SIZE AND POSITION AFTER EACH TIMESTEP
def update_moulin_wall(Mx_upstream, Mx_downstream,Mr_major, Mr_minor, dr_major,dr_minor, dGlen, dGlen_cumulative):
    Mx_upstream = Mx_upstream + dGlen - dr_major
    Mx_downstream = Mx_downstream + dGlen + dr_minor    
    if (dGlen).all == 0:
        dGlen_cumulative = 0 #this prevents from adding glen if it's been deactivated. dGlen would still be calculated in the code and unvolontary added     
    Mr_major = dGlen_cumulative - Mx_upstream #(m) relative moulin radius
    Mr_minor = Mx_downstream - dGlen_cumulative  
    #Mr_major_test = Mr_major-dr_major
    #Mr_minor_test = Mr_minor-dr_minor
    #if Mr_major != Mr_major_test:
    #    print('Houston, we have a problem in calculate_new_moulin_wall_position')
    return [Mx_upstream,Mx_downstream,Mr_major,Mr_minor]

def calculate_cumulative_dGlen(dGlen, dGlen_cumulative):

    return dGlen_cumulative + dGlen # dGlen[0] this seems to be unnecessary


#CALCULATE OUTPUTS   
def calculate_Qout(S,hw,L):
    """Discharge out of the subglacial channel"""
    return c3*S**(5/4)*np.sqrt(rhow*g*hw/L)    
        

