import numpy as np
import copy

secs_per_day = 60.0 * 60.0 * 24.0

default_params = {
    # Standard constants (normally will use default)
    "rho_i": 910.0,  # kg/m^3, ice density
    "rho_w": 1000.0,  # kg/m^3, water density
    "L_f": 3.32e5,  # J/kg, Latent heat of fusion for ice #(CT) 3.35e5 in schoof(2010)
    "g": 9.8,  # m/s^2, grav accel
    # Glen's flow law params
    "n": 3.0,  # unitless, exponent in Glen's Flow Law
    "A": 6.0e-24,  # Pa^-n /s, Constant in Glen's flow law
    # Conduit params
    "f": 0.1,  # unitless, D-W friction factor
    "C": 2.0 ** (5.0 / 4.0)
    / np.pi ** (1.0 / 4.0)
    * np.sqrt(np.pi / (np.pi + 2.0)),  # shape parameter (default is semi-circular)
    # Glacier params
    "L": 15000.0,  # m, conduit length
    "h_term": 0.0,  # m, head at terminus
    "h0_init": 500,  # m, head at moulin
    # moulin params
    "A_R": 50.0,  # m^2, xc area of moulin (reservoir)
    # Recharge params
    "R_mean": 1.0,  # m^3/s, mean discharge into moulin
    "R_func": "constant",  # type of function for recharge
    "R_period": secs_per_day,  # Period for oscillatory Recharge
    "R_min": 0.1,  # Minimum daily recharge (used for sine func)
    # Settings for 1D grid
    "nx": 50,  # Number of nodes in simulation
    "S0": 1.0,  # Conduit area in m^2
    "dt": 500,  # time in seconds
}


class one_dim_sim:
    """
    Extended 1D simulation of a subglacial conduit.
    """

    def __init__(self, params=None):
        if params is None:
            params = {}
        self.params = copy.copy(params)
        self.set_default_params()
        # Set up arrays to hold variables
        # Quantities on nodes
        self.h_node = np.zeros(self.nx)
        self.h_node[0] = self.h0_init
        self.h_node[-1] = self.h_term
        self.x = np.linspace(0, self.L, self.nx)
        self.dx = self.x[1] - self.x[0]
        # Quantities between nodes
        self.S = self.S0 * np.ones(self.nx - 1)
        self.xmid = (self.x[1:] + self.x[:-1]) / 2.0
        # Calculate ice thickness from sqrt profile and L
        H0 = 1000.0 / np.sqrt(
            30000.0
        )  # assures ice thickness of 1500 m at 60 km from edge
        self.Z = np.flip(H0 * np.sqrt(self.xmid))
        # Set up additional constants
        self.C1 = self.calc_C1()
        self.C2 = self.calc_C2()
        self.C3 = self.calc_C3()
        self.t = 0.0
        self.set_R_func()

    def set_default_params(self):
        # Define default params if they do not yet exist
        for param in default_params.keys():
            if param not in self.params:
                self.params[param] = default_params[param]
        for param in self.params:
            setattr(self, param, self.params[param])

    def set_R_func(self):
        if self.params["R_func"] == "constant":
            self.R_func = self.R_const
        if self.params["R_func"] == "sin":
            self.R_func = self.R_sin

    # Function for constant recharge
    def R_const(self, t):
        return self.params["R_mean"] * np.ones(np.size(t))

    # Function for sinusoidal recharge
    def R_sin(self, t):
        # (R_mean - R_min) sin(2 pi t/P) + R_mean
        return (self.params["R_mean"] - self.params["R_min"]) * np.sin(
            2.0 * np.pi * t / self.params["R_period"]
        ) + self.params["R_mean"]

    def calc_C1(self):
        return 1.0 / self.params["rho_i"] / self.params["L_f"]

    def calc_C2(self):
        n = self.params["n"]
        return 2.0 * self.params["A"] * n ** -n

    def calc_C3(self):
        return self.params["C"] / np.sqrt(self.params["rho_w"] * self.params["f"])

    def calc_flow(self):
        dh_tot = self.h_node[0] - self.h_node[-1]
        Q_num = self.C3 * np.sqrt(dh_tot * self.rho_w * self.g)
        Q_denom = np.sqrt(self.dx * np.sum(self.S ** (-5.0 / 2.0)))
        self.Q = Q_num / Q_denom
        # Calc heads based on Q
        self.dh = (
            self.Q ** 2
            * self.dx
            / (self.C3 ** 2 * self.S ** (5.0 / 2.0) * self.rho_w * self.g)
        )
        # Update heads
        self.h_node[1:-1] = self.h_node[0] - np.cumsum(self.dh[:-1])
        self.h_link = (self.h_node[:-1] + self.h_node[1:]) / 2.0

    def update_moulin_head(self):
        dh_dt = (self.R_func(self.t) - self.Q) / self.A_R
        self.h_node[0] += dh_dt * self.dt

    def melt_creep(self):
        rho_w = self.rho_w
        rho_i = self.rho_i
        g = self.g
        self.melt = self.C1 * self.Q * (rho_w * g * self.dh / self.dx)
        P_i = rho_i * g * self.Z
        P_w = rho_w * g * self.h_link
        self.creep = self.C2 * (P_i - P_w) ** self.n * self.S
        dS_dt = self.melt - self.creep
        dS = dS_dt * self.dt
        self.S += dS

    def run_one_step(self, dt=500.0):
        self.dt = dt
        self.calc_flow()
        self.melt_creep()
        self.update_moulin_head()
        self.t += dt
