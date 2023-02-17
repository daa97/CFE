import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.interpolate as inter
import csv
import time
from fluids import *
# import coolprop
import matplotlib as mpl
mpl.rc('font', family='Times New Roman',size="10")
mpl.rc('figure', figsize=(4.8,3.6))
mpl.rc('legend', labelspacing=0.05)
mpl.rc('savefig', dpi=800)
mpl.rc('lines', linewidth=1.2)
mpl.rc('axes', grid=True)
mpl.rc('grid', linewidth=0.25)
mpl.rc('mathtext', fontset="dejavuserif")
mpl.rc('xtick.minor', visible=True, size=1.5, width=0.5)
mpl.rc('ytick.minor', visible=True, size=1.5, width=0.5)
plt.rcParams['figure.constrained_layout.use'] =  True

opts = {
    "dim" : "Y",
    "prelim" : "y",
    "geom" : "y",
    "stations" : ["Y","y"],
    "losses" : "Y"
}

nozzle_inputs = {
        "radius ratio" : 1.1,
        "camber angle" : np.pi/6,
        "ac" : 0.25,
        "t_2c" : 0.025,
        "t_3c" : 0.012,
        "t_maxc" : 0.06,
        "dc" : 0.4,
        "sc" : 0.75,
        "setting angle" : 10/180*np.pi
    }


dynamic_turb_inputs = {
    "PR_ts" : 1.0008,
    "eta_ts" : 0.9,
    "h_0ss" : 0,
    "N_s" : 0,
    "v_s" : 0.696
}

class CFE:
    def __init__(self,static_cfe_inputs,dynamic_turb_inputs,i):

        self.i = i
        self.static_cfe_inputs = static_cfe_inputs
        """ Entry channel fluid properties """
        self.T_in = static_cfe_inputs["temp"] #Temperature [K]
        self.P_in = static_cfe_inputs["press"] * 1e6
        self.P_out = static_cfe_inputs["press"] / dynamic_turb_inputs["PR_ts"] * 1e6
        # print("CFE Outlet Pressure:",self.P_out/1e6,"[MPa]")
        self.cfe_state = H2(p=self.P_in,t = self.T_in)
        self.gamma = self.cfe_state.gamma
        T = self.T_in
        self.mu = self.cfe_state.mu
        self.rho = self.cfe_state.rho #Density, [kg/m^3]
        self.nu = self.mu/self.rho #Kinematic viscosity [m^2 s^-1]
        
        """ Entry channel geometric properties """
        self.R_i = static_cfe_inputs["inner_radius"] #CFE channel inner radius [m]
        self.R_o = static_cfe_inputs["outer_radius"] #CFE channel outer radius [m]
        self.uranium_mass = static_cfe_inputs["uranium_mass"]
        self.annulus_area = np.pi * (self.R_o**2 - self.R_i**2) #CFE channel annulus area [m^2]
        self.h = self.R_o - self.R_i #CFE channel thickness [m]
        self.eta = self.R_i/self.R_o #CFE channel radius ratio
        self.L = static_cfe_inputs["length"] #CFE Channel length [m]
        
        """ CFE Turbine requirements and inputs"""
        self.mass_flow = static_cfe_inputs["mass_flow"] #Mass Flow Rate [kg/s]
        self.omega = static_cfe_inputs["rpm"]*np.pi/30 #Angular Velocity [s^-1]
        self.work_rate = self.calc_work_rate()  #Power requirement for turbine [W]
        st_turb_inputs = self.calc_static_turb_inputs()
        # print(55 *static_cfe_inputs["rpm"]/7000) # checking the bearing resistance relative to rpm
        self.static_turb_inputs = {
            "R_i" : self.R_i,
            "R_o" : self.R_o,
            "T_01" : st_turb_inputs[0],
            "C_m1" : st_turb_inputs[1],
            "C_theta1" : st_turb_inputs[2],
            "P_01" : st_turb_inputs[3],
            "work_rate" : self.work_rate,
            "mass_flow" : self.mass_flow,
            "omega" : self.omega
        }

    def calc_mass(self):
        # TODO: add turbine mass
        PM_porosity = 0.36
        rho_SiC = 3100.2
        rho_PM = rho_SiC * (1-PM_porosity)
        m_PM = self.L * np.pi * (.049**2 - .045**2) * rho_PM

        alpha_HR1 =  (0.008/(1000-70))*(9/5)          # from NASA report, with converting F to C
        rho_case = 8070 /(1 + alpha_HR1*(self.T_in - 298.15))**3
        ri_case = self.R_i - .005                   # 5 mm thick case
        m_case = self.L * np.pi * (self.R_i**2 - ri_case**2) * rho_case

        self.mass = m_case + m_PM + self.uranium_mass
        return self.mass
    
    def calc_bearing_moment(self):
        fric_coeff = .0015
        base_load = 450
        TWR = 1.3
        diams = [.020, .020, .060]

        load = base_load + (TWR * 9.806 * self.calc_mass())
        M1 = fric_coeff * (diams[0]/2) * load       # loaded bearing
        M2 = fric_coeff * (diams[1]/2) * base_load  # unloaded bearing
        M3 = fric_coeff * (diams[2]/2) * base_load  # unloaded bearing
        self.M_bearing = M1 + M2 + M3
        return self.M_bearing

    def calc_visc_moment(self):
        Q = self.mass_flow/self.cfe_state.rho
        A = np.pi * (self.R_o**2 - self.R_i**2)
        U = Q/A
        Re_w = self.omega * self.R_i * (self.R_o - self.R_i) / self.nu
        Re_a = U * (self.R_o - self.R_i) / self.nu
        #print(f"Re w, a: {Re_w:,.1f}, {Re_a:,.1f}")
        lR = np.log10(Re_w)
        lewis_eta = 0.15999/0.22085
        lewis_G_factor = 4*np.pi*lewis_eta/((1-lewis_eta)*(1-lewis_eta**2))
        C0 = np.log10(lewis_G_factor)
        if 2600<=Re_w<=13e3:
            exp = .2005*lR**3 - 1.970*lR**2 + (7.775-1)*lR - 5.516 - C0
        elif 13e3<Re_w<=1e6:
            exp = -.006360*lR**3 + .1349*lR**2 + (.8850-1)*lR + 1.610 - C0
        else:
            raise ValueError("Rotational Reynolds number out of range!")
        lewis_Nu = 10**exp
        Nu = 1.1* lewis_Nu      # adjustment for axial flow
        #print(f"Nu: {Nu:.4f}")
        Mlam = 4*np.pi*self.cfe_state.mu*self.L*self.omega / (self.R_i**(-2) - self.R_o**(-2))
        #print(f"Wlam: {Mlam*self.omega:.4f}")
        self.M_visc = Nu*Mlam
        #print(f"W_visc: {self.M_visc*self.omega}")
        return self.M_visc

    def calc_work_rate(self):
        """Calculates the required work rate of the turbine based on viscous
        losses from shear at the CFE surface. This calculation is an 
        amalgamation of various papers of torque from taylor-couette flow.
        The correlation for non-dimensional torque comes from Lewis 1999.
        """
        M_bearings = self.calc_bearing_moment()
        M_visc = self.calc_visc_moment()        
        M = M_bearings + M_visc
        work = M * self.omega
        return work

    def calc_static_turb_inputs(self):
        U_thetaB = self.omega * self.R_i**2 / (self.R_o + self.R_i)
        U_mB = self.mass_flow / self.cfe_state.rho/self.annulus_area
        U_B = np.sqrt(U_thetaB**2 + U_mB**2)
        M_B = U_B/self.cfe_state.a
        gam = self.cfe_state.gamma
        T_01 = self.T_in * ( 1 + (gam - 1) / 2 * M_B**2 )
        P_01 = self.P_in * (T_01/self.T_in)**(gam/(gam-1))
        return [T_01,U_thetaB,U_mB,P_01]
