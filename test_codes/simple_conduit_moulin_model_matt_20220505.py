import copy
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
import matplotlib.pyplot as plt


secs_per_day = 60.*60.*24.

default_params = {
#Standard constants (normally will use default)
'rho_i':910., #kg/m^3, ice density
'rho_w':1000.,#kg/m^3, water density
'L_f':3.32e5, #J/kg, Latent heat of fusion for ice #(CT) 3.35e5 in schoof(2010)
'g':9.8, #m/s^2, grav accel

#Glen's flow law params
'n':3., #unitless, exponent in Glen's Flow Law
'A':6.e-24, #Pa^-n /s, Constant in Glen's flow law

#Conduit params
'f':0.1, #unitless, D-W friction factor
'C': 2.**(5./4.)/np.pi**(1./4.) * np.sqrt(np.pi/(np.pi + 2.)),#shape parameter (default is semi-circular)

#Glacier params
'Z_i':750., #m, ice thickness
'L': 15000., #m, conduit length
'h_term': 0., #m, head at terminus

#moulin params
'A_R': 50.,#m^2, xc area of moulin (reservoir)

#Recharge params
'R_mean':1., #m^3/s, mean discharge into moulin
'R_func':'constant', #type of function for recharge
'R_period':secs_per_day,#Period for oscillatory Recharge
'R_min':0.1,#Minimum daily recharge (used for sine func)
'overflow':True, #heads greater than overburden allowed if True
'Equil_damping_times':5,#Number of damping times to equilibrate before variable recharge starts
}

low_camp_params = {
#Glacier params
'Z_i':500., #m, ice thickness, actual is 503.
'L': 15000., #m, conduit length
'h_term': 0., #m, head at terminus
}

high_camp_params = {
#Glacier params
'Z_i':700., #m, ice thickness, actuall thickness is 712 m
'L': 25000., #m, conduit length
'h_term': 400., #m, head at terminus
}



class simple_sim():
    """
    Parameters, results, and methods of a simple conduit simulation.
    """
    def __init__(self,params=None, camp=None):
        if params is None:
            params={}
        self.params = copy.copy(params)
        #Check whether we want to run a purely nondimensional sim
        #If not, then set default params for values not given
        if 'T1' in self.params and 'T2' in self.params:
            self.nd=True
            self.T1 = self.params['T1']
            self.T2 = self.params['T2']
            if 'overflow' not in self.params:
                self.params['overflow'] = default_params['overflow']
            if 'R_mean' not in self.params:
                self.params['R_mean']=1.
            if 'R_func' not in self.params:
                self.params['R_func']= default_params['R_func']
            if 'Equil_damping_times' not in self.params:
                self.params['Equil_damping_times'] = default_params['Equil_damping_times']
            if 'R_period' not in self.params:
                self.params['R_period'] = default_params['R_period']
            #Set nondim terminal head (zero by default)
            if 'h_term_nd' in self.params:
                self.h_term_nd = self.params['h_term_nd']
            else:
                self.h_term_nd = 0.
        else:
            self.nd=False
            self.camp = camp
            self.set_default_params()
            self.calc_timescales()
            self.calc_nodim_params()
            self.h_fl =                self.params['rho_i']*self.params['Z_i']/self.params['rho_w']
            self.S_0 = (self.params['R_mean']**2.*self.params['L']/
                        (self.calc_C3()**2.*self.params['rho_i']*self.params['g']*self.params['Z_i']) )**(2./5.)
            self.h_term_nd = self.params['h_term']/self.h_fl

        #Set dummy R_func so we can evaluate timescales and equilibrium
        self.R_func = self.R_const
        self.est_damping_osc_timescales()
        self.set_R_func()
        self.calc_eq()


    def set_default_params(self):
        #Set defaults for low or high camp if specified.
        #These are set into params so that they won't be overwritten
        #from default_params.
        if self.camp == 'low':
            for param in low_camp_params:
                if param not in self.params:
                    self.params[param] = low_camp_params[param]
        if self.camp == 'high':
            for param in high_camp_params:
                if param not in self.params:
                    self.params[param] = high_camp_params[param]

        #Define default params if they do not yet exist
        for param in default_params.keys():
            if param not in self.params:
                self.params[param] = default_params[param]

    def set_R_func(self):
        if self.params['R_func'] == 'constant':
            self.R_func = self.R_const
        if self.params['R_func'] == 'sin':
            self.R_func = self.R_sin

    #Function for constant recharge
    def R_const(self, t):
        return self.params['R_mean']*np.ones(np.size(t))

    #Function for sinusoidal recharge
    def R_sin(self, t):
        equil_time = self.params['Equil_damping_times']*self.damping
        if t<equil_time:
            return self.params['R_mean']
        #(R_mean - R_min) sin(2 pi t/P) + R_mean
        if self.nd:
            #Return R_sin calculated from dimensionless parameters
            return (self.params['R_mean'] - self.params['R_min'])* np.sin(2.*np.pi*t/self.params['R_period']) + self.params['R_mean']
        else:
            #return R_sin calculated from dimensional parameters
            return (self.params['R_mean'] - self.params['R_min'])* np.sin(2.*np.pi*t*self.tau_res/self.params['R_period']) + self.params['R_mean']


    def calc_timescales(self):
        #Calculate timescales
        self.tau_res = self.calc_tau_res()
        self.tau_melt = self.calc_tau_melt()
        self.tau_creep = self.calc_tau_creep()


    def calc_nodim_params(self):
        #Make sure we've calculated timescales. If not, then calculate them.
        if not hasattr(self,'tau_res') or not hasattr(self, 'tau_melt') or not hasattr(self,'tau_creep'):
            self.calc_timescales()
        #Calculate T1 and T2 (dimensionless params)
        self.T1 = self.tau_res/self.tau_melt
        self.T2 = self.tau_res/self.tau_creep

    def est_damping_osc_timescales(self):
        y_eq = self.calc_eq()
        jac_eq = self.jac(0., y_eq)
        p = jac_eq[0,0] + jac_eq[1,1]#trace of jacobian
        q = jac_eq[0,0]*jac_eq[1,1] - jac_eq[0,1]*jac_eq[1,0]#det of jacobian
        self.alpha = p/2
        self.beta = np.imag(np.sqrt(np.complex(p**2. - 4.*q))/2.)
        self.damping = np.abs(1./self.alpha)
        if self.beta != 0: #Some cases seem to have zero beta (infinite oscillation timescale?)
            self.osc = 2.*np.pi/self.beta

    def calc_tau_melt(self):
        C1 = self.calc_C1()
        C3 = self.calc_C3()
        P_i = self.params['rho_i']*self.params['g']*self.params['Z_i']
        return (self.params['L']/P_i)**(7./5.)/C1/C3**(4./5.)/self.params['R_mean']**(1./5.)

    def calc_tau_creep(self):
        C2 = self.calc_C2()
        P_i = self.params['rho_i']*self.params['g']*self.params['Z_i']
        return 1./C2/P_i**self.params['n']#I think this is an n. should rework derivation in terms of n rather than n=3 to make sure these relations are general

    def calc_tau_res(self):
        return self.params['A_R']*self.params['rho_i']*self.params['Z_i']/self.params['rho_w']/self.params['R_mean']

    def calc_C1(self):
        return 1./self.params['rho_i']/self.params['L_f']

    def calc_C2(self):
        n= self.params['n']
        return 2.*self.params['A']*n**-n

    def calc_C3(self):
        return self.params['C']/np.sqrt(self.params['rho_w']*self.params['f'])

    def calc_eq_approx(self):
        self.h_eq_approx = 1./((self.T1/self.T2)**(5./7.) + 15./7.)
        self.S_eq_approx = self.h_eq_approx**(-2./5.)
        return (self.h_eq_approx, self.S_eq_approx)

    def calc_eq(self):
        y_eq_approx = self.calc_eq_approx()
        soln = root(self.dy_dt_no_time, y_eq_approx)
        y_eq = soln.x
        self.h_eq = y_eq[0]
        self.S_eq = y_eq[1]
        return y_eq

    #This swaps the order of t and y for compatibility with root()
    def dy_dt_no_time(self,y):
        return self.dy_dt(0,y)#set t to zero

    def dy_dt(self, t, y):
        h = y[0]
        S = y[1]
        if not self.params['overflow']:
            if h>1:
                h=0.999
        dh_dt = self.R_func(t)/self.params['R_mean'] - S**(5./4.)*h**.5
        dS_dt = self.T1*S**(5./4.)*np.abs((h-self.h_term_nd))**(3./2.) - self.T2*S*(1.-h)**3.
        if not self.params['overflow']:
            if h>0.990:
                if dh_dt>0:
                    dh_dt=0.
        return (dh_dt, dS_dt)

    #Function to calculate Jacobian
    def jac(self, t, y):
        h = y[0]
        S = y[1]
        T1 = self.T1
        T2 = self.T2
        j11 = -0.5*h**-0.5*S**(5./4.)
        j12 = -(5./4.)*h**0.5*S**0.25
        j21 = (3./2.)*T1*S**(5./4.)*h**0.5 + 3.*T2*(1.-h)**2.*S
        j22 = (5./4.)*T1*S**0.25*h**1.5 - T2*(1.-h)**3.
        return np.array([[j11,j12],[j21,j22]])

    def run(self, run_params={}):
        if not 't_i' in run_params:
            run_params['t_i']= 0
        if not 't_f' in run_params:
            run_params['t_f'] = 30
            if self.params['R_func'] == 'sin':
                run_params['t_f'] += self.damping*self.params['Equil_damping_times']
        if not 'S_0' in run_params:
            run_params['S_0'] = self.S_eq*0.5
        if not 'h_0' in run_params:
            run_params['h_0'] = self.h_eq*0.5
        if not 'method' in run_params:
            run_params['method'] = 'LSODA'
        if not 'atol' in run_params:
            run_params['atol'] = 1e-6
        if not 'rtol' in run_params:
            run_params['rtol'] = 1e-3
        if not 'max_step' in run_params:
            run_params['max_step'] = 0.1
        sol = solve_ivp(self.dy_dt,
            (run_params['t_i'], run_params['t_f']),
            (run_params['h_0'],run_params['S_0']),
            jac=self.jac,
            method=run_params['method'],
            atol=run_params['atol'],
            rtol = run_params['rtol'],
            max_step = run_params['max_step'])
        self.sol = sol
        self.h_nd = sol.y[0]
        self.S_nd = sol.y[1]
        self.t_nd = sol.t
        self.Q_nd = self.S_nd**(5./4.)*self.h_nd**0.5
        if not self.nd:
            self.h = self.h_nd*self.h_fl
            self.S = self.S_nd*self.S_0
            self.t = self.t_nd*self.tau_res
            self.Q = self.Q_nd*self.params['R_mean']
        return sol

    def plot_sim(self, nd=True, separate=False, filename=None, close_plots=False, tmin=None, as_overburden=False):
        plt.figure()
        if nd:
            h = self.h_nd
            S = self.S_nd
            t = self.t_nd
        else:
            h = self.h
            if as_overburden:
                h = self.h_nd
            S = self.S
            t = self.t/secs_per_day
        if not separate:
            plt.plot(t,h,t,S)
            if tmin != None:
                plt.xlim([tmin, max(t)])
            plt.legend(['h','S'])
            if nd:
                plt.xlabel('$t$')
                plt.ylabel('$h$, $S$')
            else:
                plt.xlabel(r'$t\,{\rm (days)}$')
                plt.ylabel(r'$h\,{\rm(m)},\, S\,{\rm (m^2)}$')
            if filename != None:
                plt.savefig(filename)
            if close_plots:
                plt.close()

        else:
            plt.plot(t,h)
            if tmin != None:
                plt.xlim([tmin, max(t)])
            if nd:
                plt.xlabel('t')
                plt.ylabel('h')
            else:
                plt.xlabel(r'$t\,{\rm (days)}$')
                plt.ylabel(r'$h\,{\rm(m)}$')
            if filename != None:
                plt.savefig(filename[:-4]+'-head'+filename[-4:])
            if close_plots:
                plt.close()
            plt.figure()
            plt.plot(t,S)
            if tmin != None:
                plt.xlim([tmin, max(t)])
            if nd:
                plt.xlabel('t')
                plt.ylabel('S')
            else:
                plt.xlabel(r'$t\,{\rm (days)}$')
                plt.ylabel(r'$S\,{\rm(m^2)}$')
            if filename != None:
                plt.savefig(filename[:-4]+'-area'+filename[-4:])
            if close_plots:
                plt.close()


    def plot_phase(self, nd=True, filename=None, close_plots=False):
        plt.figure()
        if nd:
            h = self.h_nd
            S = self.S_nd
        else:
            h = self.h
            S = self.S
        plt.plot(h,S)
        if nd:
            plt.xlabel('$h$')
            plt.ylabel('$S$')
        else:
            plt.xlabel(r'$h\, {\rm (m)}$')
            plt.ylabel(r'$S\, {\rm (m^2)}$')
        if filename != None:
            plt.savefig(filename)
        if close_plots:
            plt.close()
