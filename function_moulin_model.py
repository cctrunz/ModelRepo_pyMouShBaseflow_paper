'''
Moulin physical model Functions for run_moulin_model.

Code developped and writen in Matlab by Kristin Poinar and Lauren Andrews, 2019-2020.
Translation in Python by Celia Trunz, 2020.
'''

import numpy as np


'''Default constants'''
T0 = 273.15 #Kelvin
rhoi = 910 #kg/m3; Ice density
rhow = 1000 #kg/m3; Water density
ki = 2.1 #J/mKs
cp = 2115 #J/kgK
Lf = 335000 #J/kg; Latent heat of fusion
g = 9.8 #m/s2; Gravity
E = 5e9 #Pa; Young's elastic modulus (Vaughan 1995)
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

def initiate_moulin_radius(z,r_minor_inital,r_major_inital):  
    '''
    - r_minor is the radius that controls the circle and the small axis of the elipse
    and r_major control the large axis of the elipse
    '''
    #Define initial moulin characteristics
    # create z with generate_grid_z()
    Mr_minor = r_minor_inital * np.ones(len(z))
    Mr_major = r_major_inital * np.ones(len(z))
    return [Mr_minor, Mr_major]
    
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
    #dz:  vertical spacing between the vertical point of calculation in the moulin
    return np.arange(0,H+1,dz) #(m) vertical profile fo moulin point of calculation 
        #!! we add +1 at xmax because python consider it not inclusive. 
        #This way, entered parameters are going to match matlab input

def generate_time(dt,tmax_in_day):
    tmax_in_second = tmax_in_day*24*3600
    time_vector = np.arange(dt,tmax_in_second+1,dt)#+1 is important to match matlab
    return [time_vector, tmax_in_second]
        
#def generate_vector_time(number_of_days=20,timestep_in_seconds=300):
    
       
def set_ice_temperature(x,z, type='Temperate'):
    '''Define ice temperature
    - Create x with generate_grid_x
    - Create z with generate_grid_z
    '''   

    #define temperature in the glacier, far from the moulin
    if type == 'Temperate':
        T_far = T0 * np.ones(len(z))
    if type == 'Cool':
        Tmin = -5;
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

def flow_law_parameter(T,Pi_z):
    '''Glen's Flow Law -- Cuffey and Paterson Eqn. 3.35 
    -(CT)Verify that this is true. Found info in matlab code'''   
    Tfrac = 1/(T + a*Pi_z) - 1/(Tstar + a*Pi_z)
    Qc = Qcless * np.ones(len(T))
    #Qc[T>Tstar] = Qcmore
    return Astar * np.exp(-Qc/R * Tfrac)
    
def ice_pressure_at_depth(H,z): #ask Kristin how to call that
    ''' Ice pressure in function of depth
    defined as P in matlab code'''
    return rhoi*g*(H-z)

def water_pressure_at_depth(hw,z):
    return rhow*g*(hw-z)
  
def locate_water(hw,z):
    '''Boolean index of position of z that are underwater. 
    in the Matlab code, it correspond to "wet" '''
    return z <= hw
    
#def stress_cryo(Pi_z):
#    '''Ice hydrostatic stress (INWARD: Negative)'''
#    return - Pi_z
    
#def stress_hydro(hw,z):
#    '''Water hydrostatic stress (OUTWARD: Positive)'''
#    stress_hydro = rhow*g*(hw-z) #Calculate water pressure at a certain depth in the moulin
#    return stress_hydro(np.invert(locate_water()))

def S_moulin_at_h(h, Mr_minor, Mr_major):
    '''Calculate cross-section area of the moulin at the water level. 
    The cross-section area is composed of two part: 1 demi-circle and 1 demi-elipse
    '''
    #cross-section area of demi-circle. 
    S_demi_circle = (Mr_minor**2 * np.pi)/2
    S_demi_elipse = (Mr_minor * Mr_major * np.pi)/2
    return S_demi_circle+S_demi_elipse
    
#def Qin_func(type='Sinusoidal_Celia'): #there should be a better way to do this...
#    if type=='Constant':
#        return Qin_constant
#    ##R_sin
#    if type=='Sinusoidal_Celia':
#        return Qin_Sinusoidal_Celia()
    
def Qin_constant(Qin=3):
    return Qin

def Qin_Sinusoidal_Celia(t,Qin_mean=3, Qin_min=2, period=24*3600):

    return (Qin_mean- Qin_min) * np.sin(2*np.pi*t/period) + Qin_mean
#MAIN FUNCTIONS
    

#CREEP
def creep_moulin(Mr,dt,T,Pi_z,stress):
    ''' Creep closure of a water-filled borehole   
    Based on boreholeclosure/HomeworkProblem_Vostok3G.m which Krisin Poinar did in 2013 for crevasse model    
    Borehole 3G at Vostok by Blinov and Dmitriev (1987) and Salamatin et al (1998) from Table 4 in Talalay and Hooke, 2007 (Annals)    
    .. code writen in Matlab by Kristin Poinar
    '''  
    T_mean = np.mean(T, axis=1) # mean of each row of tempearture in the x axis in matlab, it's written mean(T,2)
    A = flow_law_parameter(T_mean, Pi_z) #
    #total stress
    sigma_z = stress['cryo'] + stress['hydro'] 
    
    epsilon_dot = E*A*(sigma_z/3)**3 #boreholeclosure/HomeworkProblem_Vostok3G.m divided A by 5 in order to match measured Antarctic BH closure rates
    #Creep closure rate
    return Mr*np.exp(epsilon_dot*dt)-Mr

def function_subglacial_schoof(t,y,hw,Sm_h,Pi,L,Qin):
    '''Define dimensionalized function for subglacial channel and water level in moulin
    input:
    - t = time vector
    - matrice containing initial head and channel cross-section area
    - L = conduit length
    - Pi = ice pressure
    output:
    - head timeserie
    - channel cross-section area time serie  
    .. code from Celia Trunz, from Schoof 2010, subglacial channel only, without the cavities'''

    # print(f'y: {y}')
    h = y[0] #(–) S*: Non dimentionalized channel cross-section area
    S = y[1] #(–) h*: Non dimentionalized moulin head

    #Head partial differential equation
    dhdt = 1/Sm_h * ( Qin - c3*S**(5/4)*np.sqrt((rhow*g*hw)/L) )
    #Channel cross section area partial differential equation
    dSdt = c1 * c3 * S**(5/4) * ((rhow*g*h)/L)**(3/2) - c2 * ( Pi - rhow*g*h )**n * S
    return [dhdt, dSdt]

