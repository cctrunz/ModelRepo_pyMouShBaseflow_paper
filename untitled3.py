import numpy as np
import function_moulin_model as fmm
from scipy.integrate import solve_ivp
from functools import partial


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


'''Default parameter'''
#R0 = 2 #(m) Initial moulin radius
H = 500 #(m) Ice thickness
L = 10000 #(m) Subglacial channel length
tmax_in_day = 20 #(days) Maximum time to run
dt = 300 #(s) timestep
r_minor_initial = 2 #(m)
r_major_initial = 2 #(m)
xmax    = 30 # 80 #(m) how far away from moulin to use as infinity
hw = H/2 #(m)Initial water level
S = 2 #(m)
Qin = 3 #m3/s


def function_subglacial_schoof(t,y,Ms,z,Pi,L,Qin):
    '''Define dimensionalized function for subglacial channel and water level in moulin
    input:
    - t = time vector
    - matrice containing initial head and channel cross-section area
    - Ms = moulin cross-section area at each depth [array]
    - L = conduit length
    - Pi = ice pressure
    output:
    - head timeserie
    - channel cross-section area time serie  
    .. code from Celia Trunz, from Schoof 2010, subglacial channel only, without the cavities'''
    hw = y[0] #(m) moulin head initial value
    S = y[1] #(m) channel cross-section area initial value
    
    Ms_hw = np.interp(hw,z,Ms)#find Ms value by interpolating Ms value at the level of the water
    dhwdt = 1/Ms_hw* ( Qin - c3*S**(5/4)*np.sqrt((rhow*g*hw)/L) )# Moulin head ODE
    dSdt = c1 * c3 * S**(5/4) * ((rhow*g*hw)/L)**(3/2) - c2 * ( Pi - rhow*g*hw )**n * S# Channel cross-section area ODE
    return [dhwdt, dSdt]

[z, nz, dz] = fmm.generate_grid_z(H) #default is 1m spacing # should this mimic the x grid??
Pi = fmm.ice_pressure_at_depth(H,0) #total ice pressure
[Mr_minor, Mr_major] = fmm.initiate_moulin_radius(z,r_minor_initial,r_major_initial) #
Ms = (np.pi*Mr_minor**2)/2 + (np.pi*Mr_minor*Mr_major)/2 #Moulin cross-section area 
[time,tmax_in_second] = fmm.generate_time(dt,tmax_in_day)

'''Calculate water level'''
function_compressed = partial(function_subglacial_schoof, Ms, z, Pi, L, Qin) 
sol = solve_ivp(#function_compressed,
                function_subglacial_schoof,
                #fun=lambda t,y: fmm.function_subglacial_schoof(t,y, Ms, z, Pi, L, Qin),
                [0, dt], #initial time and end time. We solve for one timestep.
                [hw,S], #initial head and channel cross-section area. Uses values in previous timestep.
                args=(Ms,z,Pi,L,Qin),
                method = 'LSODA' #solver method
                # atol = 1e-6, #tolerance. can make it faster
                # rtol = 1e-3,
                #max_step = 10
                )

hw = sol.y[0][-1]  #(m) moulin water level
S = sol.y[1][-1] #(m) Channel cross-section