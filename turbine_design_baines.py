import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.interpolate as inter
import csv
from fluids import *

"""
Syntax for thermo calculations

from fluids import *


H2 = fluid("Hydrogen", prop_files)

state = H2(h=H_1, P=P_1)
state = H2(s=s_1,P=P_1)

T = state.T
or 
h_1 = state.h #J/kg
Cp_1 = state.cp 

"""
class CFE:
    def __init__(self,static_cfe_inputs,dynamic_turb_inputs,i):
        self.i = i
        self.static_cfe_inputs = static_cfe_inputs
        """ Entry channel fluid properties """
        self.T_in = static_cfe_inputs["temp"] #Temperature [K]
        self.P_in = static_cfe_inputs["press"] * 1e6
        self.P_out = static_cfe_inputs["press"] / dynamic_turb_inputs["PR_ts"] * 1e6
        # print("CFE Outlet Pressure:",self.P_out/1e6,"[MPa]")
        self.cfe_state = H2(P=self.P_in,T = self.T_in)
        self.gamma = self.cfe_state.gamma
        T = self.T_in
        self.mu = eval("-0.00000000000144 *T**2 + 0.0000000169 *T+ 0.00000464")
        self.rho = self.cfe_state.rho #Density, [kg/m^3]
        self.nu = self.mu/self.rho #Kinematic viscosity [m^2 s^-1]
        
        """ Entry channel geometric properties """
        self.R_i = static_cfe_inputs["inner_radius"] #CFE channel inner radius [m]
        self.R_o = static_cfe_inputs["outer_radius"] #CFE channel outer radius [m]
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
            "T_01" : st_turb_inputs[0],
            "C_m1" : st_turb_inputs[1],
            "C_theta1" : st_turb_inputs[2],
            "P_01" : st_turb_inputs[3]
        }
        
    def calc_work_rate(self):
        """Calculates the required work rate of the turbine based on viscous
        losses from shear at the CFE surface. This calculation is an 
        amalgamation of various papers of torque from taylor-couette flow.
        The correlation for non-dimensional torque comes from Lewis 1999.
        
        """
        
        Re = self.omega * self.R_i * (self.R_o - self.R_i) / self.nu #Inner cylinder Reynolds number
        
        if Re <= 13000 and Re > 2600:
            G = 10**(0.2005*(np.log10(Re))**3 -1.970 * (np.log10(Re))**2 + 7.775 * (np.log10(Re)) - 5.516 )
            
        else: 
            G = 10**(-0.006360 * (np.log10(Re))**3 + 0.1349 * (np.log10(Re))**2 + 0.8850 * (np.log10(Re)) + 1.610 )

        eta =15.999/22.085 
    
        G_lam = 4*np.pi * eta /( (1-eta)*(1-eta**2))*Re
        Nu_omega = G/G_lam
        G_lam = 2 * self.eta /( (1-self.eta)*(1-self.eta**2))*Re
        G = Nu_omega * G_lam
        T_1 = G * 2*np.pi*self.L*self.rho*self.nu**2
        T_lam = 4*np.pi*self.L*self.mu*self.omega/(self.R_i**-2 - self.R_o**-2) * 1.1 # add 10% to compensate for axial throughflow (Yamada, 1962)
        T_2 = T_lam * Nu_omega 
        T_bearings = 0.04358825
        # print("Viscous torque:",T_1)
        T_1 += T_bearings
        work = T_1*self.omega
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

def eval_scalar_input(inpt):
    if isinstance(inpt, str):
        return eval(inpt)
    else:
        return inpt 

class station:
    def __init__(self,turb):
        pass
class whitfield_turbine:
    def __init__(self,CFE,static_turb_inputs,dynamic_turb_inputs,i):
        self.CFE = CFE
        self.i = i
        self.state_01 = H2(t=static_turb_inputs["T_01"], p = self.CFE.P_in)
        self.T_01 = static_turb_inputs["T_01"]
        self.h_01 = self.state_01.h
        self.delta_W = self.CFE.work_rate / self.CFE.mass_flow
        self.S = self.delta_W / self.h_01 #Power ratio
        self.T_03 = (1 - self.S) * self.T_01

class turbine:
    def __init__(self,CFE,static_turb_inputs,dynamic_turb_inputs,i):
        self.CFE = CFE
        self.static_turb_inputs = static_turb_inputs
        self.i = i
        print(f'Turbine iteration {self.i}')
        #self.number_stations = static_turb_inputs["number of stations"]
        self.stations = self.make_stations()
        
        #self.eta_tt_guess = dynamic_turb_inputs["total to total efficiency guess"]
        self.eta_ts_guess = dynamic_turb_inputs["eta_ts"]
        self.PR_ts = dynamic_turb_inputs["PR_ts"]

        self.deltah_0 = self.CFE.work_rate/self.CFE.mass_flow
        self.T_01 = static_turb_inputs["T_01"]
        self.state_01 = H2(t=self.T_01, p = static_turb_inputs["P_01"])
        self.h_01 = self.state_01.h
        self.h_05 = self.state_01.h - self.deltah_0
        self.state_05 = H2(h = self.h_05, P=static_turb_inputs["P_01"]/self.PR_ts)
        self.epsilon_b = 0.4e-3
        self.Q = self.CFE.mass_flow / self.state_05.rho
        # print(self.state_01)
        # print(self.state_05)
        if self.i == 1:
            self.eta_it_output = self.iterate_eta_ts()
            self.eta_ts = self.eta_it_output[0]
            self.h_0ss = self.eta_it_output[1]
            self.N_s = self.eta_it_output[2]
        else:
            self.eta_ts = dynamic_turb_inputs["eta_ts"]
            self.PR_ts = dynamic_turb_inputs["PR_ts"]
            self.h_0ss = dynamic_turb_inputs["h_0ss"]
            self.N_s = self.CFE.omega * np.sqrt(self.Q) / (self.h_0ss**(3/4))

        """Preliminary Calcs: Rotor"""

        """Isentropic Values"""
        self.v_s = 0.6913 # 0.737 * self.N_s**0.2
        self.C_0 = np.sqrt(2 * self.h_0ss)
        """Total and static states at station 4"""
        self.U_4 = self.C_0 * self.v_s # Blade spped [m/s]
        self.r_4 = self.U_4/self.CFE.omega # Rotor inlet radius [m]
        self.alpha_4 = 90 - 1* (10.8 + 14.2*self.N_s**2) # Rotor inlet angle [degrees]
        self.alpha_4_rad = self.alpha_4 / 180 * np.pi # Rotor inlet angle [radians]
        self.P_04 = self.CFE.P_in - self.state_01.rho * self.h_0ss*(1-self.eta_ts)/4 # Rotor inlet total pressure [Pa]
        self.C_theta4 = self.U_4 * self.eta_ts / (2*self.v_s**2) # Rotor inlet tangiential velocity [m/s]
        self.C_m4 = self.C_theta4 / np.tan(self.alpha_4*np.pi/180)
        self.C_4 = np.sqrt(self.C_m4**2 + self.C_theta4**2)
        self.W_theta4 = self.C_theta4-self.U_4
        self.W_m4 = self.C_m4
        self.beta_4_rad = np.arctan(self.W_theta4/self.W_m4) 
        self.beta_4 = self.beta_4_rad * 180 / np.pi
        self.W_4 = np.sqrt(self.C_m4**2 + (self.C_theta4-self.U_4)**2)
        self.h_04 = self.state_01.h
        self.state_04 = H2(h=self.h_04,p=self.P_04)
        self.s_4 = self.state_04.s
        self.h_4 = self.h_04 - 0.5 * self.C_4**2  
        self.state_4 = H2(s=self.s_4,h=self.h_4)
        self.b_4 = self.CFE.mass_flow / (2*np.pi*self.r_4*self.state_4.rho*self.C_m4)
        self.A_4 = 2 * np.pi * self.r_4 *self.b_4
        self.M_4 = self.C_4 / self.state_4.a
        self.M_4_rel = self.W_4 / self.state_4.a

        """Total and static states at station 5"""
        self.r_h5 = 0.185*self.r_4 #Rotor outlet hub radius
        self.zeta = 1 + 5 * (self.b_4/self.r_4)**2 #Rotor meridional velocity ratio
        self.C_m5 = self.zeta*self.C_m4 #Rotor outlet meridional velocity
        self.C_5 = self.C_m5 #Rotor absolute velocity
        self.h_5 = self.h_05 - 0.5*self.C_5**2
        self.s_5 = self.state_05.s
        self.state_5 = H2(h=self.h_5,s=self.s_5)
        self.A_5 = self.CFE.mass_flow/(self.state_5.rho * self.C_m5)
        self.r_s5 = np.sqrt((self.A_5/np.pi)+self.r_h5**2) #Rotor outlet shroud radius
        self.r_5 = (self.r_s5+self.r_h5)/2 #Rotor outlet mean radius
        self.b_5 = self.r_s5-self.r_h5 #Rotor outlet blade height
        self.U_5 = self.CFE.omega * self.r_5
        self.W_m5 = self.C_m5
        self.W_theta5 = - self.U_5
        self.W_5 = np.sqrt(self.W_m5**2 + self.W_theta5**2)
        self.beta_5_rad = np.arctan(self.W_theta5/self.W_m5)
        self.beta_5 = self.beta_5_rad * 180 / np.pi
        self.P_05 = self.state_05.p
        self.state_5ss = H2(s = self.state_01.s,p = self.state_5.p)
        self.M_5 = self.C_5/ self.state_5.a
        self.M_5_rel = self.W_5 / self.state_5.a
        self.beta_s5 = np.arctan(self.r_5/self.r_s5 * self.beta_5_rad) * 180 / np.pi
        self.beta_h5 = np.arctan(self.r_5/self.r_h5 * self.beta_5_rad) * 180 / np.pi

        self.z_r = 1.5 * (self.b_5)
        self.n_r = np.round(np.pi/30*(110-self.alpha_4)*np.tan(self.alpha_4*np.pi/180)) #Glassman 1972
        self.q_5 = 2 * np.pi * self.r_5 / self.n_r
        self.o_5 = self.q_5 * self.C_m5 / self.W_5
        self.r_3 = self.r_4 + 2 * self.b_4 * np.cos(self.alpha_4*np.pi/180)
        self.r_2 = 1.1 * self.r_3
        self.t_lead = 0.04 * self.r_4
        self.t_trail = 0.02 * self.r_4

        """Preliminary Calcs: Nozzle"""
        self.N_n = self.n_r + 2 # Number of nozzle blades guess
        self.C_theta3 = self.r_4 / self.r_3
        self.h_03 = self.state_01.h
        self.P_03 = self.P_04
        self.b_3 = self.b_4
        self.state_03 = H2(h=self.h_03,p=self.P_03)
        rho_3_guess = 0
        rho_3_check = self.state_03.rho
        while rho_3_check - rho_3_guess != 0:
            rho_3_guess = rho_3_check
            C_m3 = self.CFE.mass_flow / (2 * np.pi * self.r_3 * self.b_3 * rho_3_guess)
            C_3 = np.sqrt(C_m3**2 + self.C_theta3**2)
            h_3 = self.h_03 - 0.5 * C_3**2
            s_3 = self.state_03.s
            state_3 = H2(s=s_3,h=h_3)
            rho_3_check = state_3.rho

        self.C_m3 = C_m3
        self.state_3 = state_3
        self.C_3 = C_3
        self.alpha_3_rad = np.arctan(self.C_theta3/self.C_m3)
        self.alpha_3 = self.alpha_3_rad * 180 / np.pi
        self.d_3 = 2 * np.pi * self.r_3 / self.N_n # Pitch is usually s but since s is entropy I am using d
        self.o_3 = self.d_3 * np.sin(self.alpha_3_rad)

        """Preliminary Calcs: The Bowl"""
        self.r_1o = self.CFE.R_o
        self.r_1i = self.r_1o - self.b_3
        self.r_cbo = self.r_1o - self.r_3
        self.r_cbi = self.r_1i - self.r_3
        # print(self.r_cbo)
        # print(self.r_cbi)

        """Efficiency Correction Calulations"""
        self.h_loss = self.calc_losses()
        self.eta_ts_loss = self.deltah_0/(self.deltah_0 + self.h_loss)
        self.flow_coef = self.C_m5 / self.U_4
        self.load_coef = self.deltah_0 / (self.U_4**2)
        self.D_s = 2 * self.r_4 * (self.h_0ss)**0.25 / self.Q**0.5

    def nozzle_blade_design(self,s_div_c = 0.75, theta = 0, a_div_c = 0.5, t_2_div_c = 0.025,
        t_3_div_c = 0.012, t_max_div_c = 0.06, d_div_c = 0.4,radius_ratio=1.1,camber_angle=0):
        N_p_n = 30
        x_div_c = np.linspace(0,1,N_p_n)
        y_div_c_guess = 1
        y_div_c_check = 0
        b_div_c = (np.sqrt(1 + (4 * np.tan(camber_angle))**2 * (a_div_c * - a_div_c**2 - 3/16))-1) / (4 * np.tan(camber_angle))
        for x in x_div_c:
            while y_div_c_guess != y_div_c_check:
                y_div_c_guess = y_div_c_check
                y_div_c_num = (x * (1 - x))
                y_div_c_den_1 = (1 - 2 * a_div_c)**2 / (4 * b_div_c**2) * y_div_c_guess
                y_div_c_den_2 = (1 - 2 * a_div_c) / (b_div_c) * x
                y_div_c_den_3 = (1 - 4 * a_div_c) / (4 * b_div_c)
                y_div_c_check = y_div_c_num / (y_div_c_den_1 + y_div_c_den_2 - y_div_c_den_3)
        
    def turbine_feasibility_checks(self):
        pass

    def make_stations(self):#use this for making an h-s diagram later
        """Could I technically make all variables equal to -1 and then
        calculate them as their dependents become available?"""
        pass    

    def iterate_eta_ts(self):
        rho_05 = self.state_05.rho
        Q = self.CFE.mass_flow / rho_05
        eta_ts = 0
        """Eta Iteration"""
        while round(eta_ts,7) != round(self.eta_ts_guess,7):
            if eta_ts != 0:
                self.eta_ts_guess = eta_ts

            h_0ss = self.deltah_0 / self.eta_ts_guess
            
            N_s = self.CFE.omega * np.sqrt(Q) / (h_0ss**(3/4))

            eta_ts = 0.87 - 1.07 * (N_s - 0.55)**2 - 0.5 * (N_s - 0.55)**3
            # print("Calculated eta:", eta_ts)
            # print("Previous eta",self.eta_ts_guess)

        h_0ss = self.deltah_0 / eta_ts
        N_s = self.CFE.omega * np.sqrt(Q) / (h_0ss**(3/4))
        return [eta_ts,h_0ss,N_s]
    
    def calc_losses(self):
        """Values required for loss calculations"""
        beta = np.arctan(0.5 * (np.tan(self.beta_4_rad) + np.tan(self.beta_5_rad)))
        chord_len = self.z_r / np.cos(beta)
        r_t = self.r_5 #estimate the throat radius as the mean rotor outlet radius
        o_t = self.o_5
        b_t = o_t
        beta_t = 0.8 * self.beta_5_rad
        W_t = 0.7*self.W_5
        L_h = np.pi / 4 *((self.z_r - (self.b_4/2)) + (self.r_4 - r_t - b_t/2))
        D_h_4 = (4 * np.pi * self.r_4 * self.b_4) / (2 * np.pi * self.r_4 + self.n_r * self.b_4 )
        D_h_5 = (2 * np.pi * (self.r_s5**2 - self.r_h5**2)) / ( np.pi * (self.r_s5 - self.r_h5) + self.n_r * self.b_5)
        D_h = 0.5 * (D_h_4 + D_h_5)

        pass_val_check = (self.r_4 - r_t) / b_t
        if pass_val_check >= 0.2:
            m_f = 1
        elif pass_val_check < 0.2:
            m_f = 2
        else: 
            m_f = 0
            ValueError("Passage losses check value is incorrect")
        
        rho_avg = (self.state_5.rho + self.state_4.rho)/2
        mu_avg = (self.state_5.mu + self.state_4.mu)/2
        c_avg = (self.C_4 + self.C_5)/2
        Re_4_avg = rho_avg * c_avg * self.r_4 / mu_avg

        if Re_4_avg >= 10**5:
            K_f = 0.0102 * (self.epsilon_b/self.r_4)**0.1 / Re_4_avg**0.2
        else:
            K_f = 3.7 * (self.epsilon_b/self.r_4)**0.1 / Re_4_avg**0.5

        M_5_rel = self.M_5_rel

        """Incidence Losses"""
        L_i_a = (0.75 / (2 * (1 - self.M_4**2)))
        L_i_b = ((self.C_m4 / self.C_m5) * (self.C_m4 / self.U_4) * np.tan(self.alpha_4*np.pi/180) - 1)**2
        L_i =  L_i_a * L_i_b * self.U_4**2
        # print("Incidence Losses:",L_i)

        """Passage Losses"""
        L_p_f = L_h/D_h
        L_p_sf = 0.68 * (1 - (r_t/self.r_4)**2) * (np.cos(beta_t) / (b_t/chord_len))
        L_p = m_f * 0.11 * (L_p_f + L_p_sf) * (self.W_4**2 + W_t**2)/2
        # print("Passage Losses:",L_p)

        """Windage Losses"""
        L_w = K_f * (rho_avg * self.U_4**3 * self.r_4**2) / (2 * self.CFE.mass_flow * self.W_5**2)
        # print("Windage Losses:",L_w)

        """Trailing Edge Losses"""
        g = 9.81
        dP_part2 = (self.n_r * self.t_lead / (np.pi * (self.r_h5 + self.r_s5) * np.cos(self.beta_5_rad)))**2
        deltaP_0_rel = self.state_5.rho * self.W_5**2 / (2 * g) * dP_part2
        
        L_t_den = self.state_5.p * (1 + (self.W_5**2/(2 * self.state_5.t * self.state_5.cp)))**(self.state_5.gamma/(self.state_5.gamma-1))
        L_t = 2 /(self.state_5.gamma * M_5_rel**2 )  * (deltaP_0_rel/L_t_den)
        # print("Trailing edge losses:",L_t)

        """Exit Energy Loss"""
        L_e = self.C_5**2 / 2
        # print("Exit Energy Losses:",L_e)

        """Nozzle Losses"""
        L_n = self.h_0ss * (1 - 0.975)

        h_losses = L_e + L_i  + L_p + L_w + L_t + L_n
        return h_losses

    def print_turbine(self,opts):
        turb_title = f'Turbine iteration: {self.i}\n\n'
        prelim_title = f'Preliminary Calculated Values:\n\n'
        dim_title = f'Non-dimensional Parameters:\n\n'
        geom_title = f'Turbine Geometry:\n\n'
        stn1_title = f'Station 1:\n\n'
        stn2_title = f'Station 2:\n\n'
        stn3_title = f'Station 3:\n\n'
        stn4_title = f'Station 4:\n\n'
        stn5_title = f'Station 5:\n\n'
        stn6_title = f'Station 6:\n\n'

        dim_str = f' \
Total-to-static efficicency: {self.eta_ts}\n \
Stage loading: {self.load_coef}\n \
Flow coefficient: {self.flow_coef}\n \
Specific speed: {self.N_s}\n \
Specific diameter: {self.D_s}\n \
Turbine velocity ratio: {self.v_s}\n \
Total-to-static pressure ratio: {self.PR_ts}\n\n'

        prelim_str = f' \
Isentropic enthalpy drop: {self.h_0ss} [kJ kg^-1]\n \
Total-to-total enthalpy drop: {self.deltah_0} [kJ kg^-1]\n \
Spouting Velocity: {self.C_0} [m s^-1]\n \
Turbine rotor inlet blade speed: {self.U_4} [m s^-1]\n\n'

        geom_str = f' \
Rotor inlet radius: {self.r_4*1000} [mm]\n \
Rotor inlet blade height: {self.b_4*1000} [mm]\n \
Rotor inlet blade angle: {self.beta_4} [degrees]\n \
Rotor outlet hub radius: {self.r_h5 * 1000} [mm]\n \
Rotor outlet shroud radius: {self.r_s5 * 1000} [mm]\n \
Rotor outlet blade height: {self.b_5 * 1000} [mm]\n \
Rotor outlet blade angle: {self.beta_5} [degrees]\n \
Rotor outlet hub blade angle: {self.beta_h5} [degrees]\n \
Rotor outlet shroud blade angle: {self.beta_s5} [degrees]\n \
Rotor axial length: {self.z_r * 1000} [mm]\n \
Rotor leading edge thickness: {self.t_lead * 1000} [mm]\n \
Rotor trailing edge thickness: {self.t_trail * 1000} [mm]\n \
Number of rotor blades: {self.n_r}\n \
Stator outlet radius: {self.r_3 * 1000} [mm]\n \
Stator inlet radius: {self.r_2 * 1e3} [mm]\n\n'

        stn4_str = f' \
Rotor inlet stagnation pressure: {self.P_04/1000} [kPa]\n \
Rotor inlet radius: {self.r_4*1000} [mm]\n \
Absolute meridional velocity: {self.C_m4} [m s^-1]\n \
Absolute tangential velocity: {self.C_theta4} [m s^-1]\n \
Absolute velocity: {self.C_4} [m s^-1]\n \
Absolute Mach no.: {self.M_4} [m s^-1]\n \
Relative tangential velocity: {self.W_theta4} [m s^-1]\n \
Relative velocity: {self.W_4} [m s^-1]\n \
Relative Mach no.: {self.M_4_rel} [m s^-1]\n \
Rotor inlet absolute flow angle: {self.alpha_4} [degrees]\n \
Rotor inlet relative flow angle: {self.beta_4} [degrees]\n\n'

        stn5_str = f' \
Rotor outlet stagnation pressure: {self.P_05/1000} [kPa]\n \
Rotor outlet radius: {self.r_5*1000} [mm]\n \
Absolute meridional velocity: {self.C_m5} [m s^-1]\n \
Absolute tangential velocity: {0} [m s^-1]\n \
Absolute velocity: {self.C_5} [m s^-1]\n \
Absolute Mach no.: {self.M_5} [m s^-1]\n \
Relative Mach no.: {self.M_5_rel} [m s^-1]\n \
Rotor outlet absolute flow angle: {0} [degrees]\n \
Rotor outlet relative flow angle: {self.beta_5} [degrees]\n\n'

        strings = {
            "dim" : dim_title + dim_str,
            "prelim" : prelim_title + prelim_str,
            "geom" : geom_title + geom_str,
            "stations" : [stn4_title + stn4_str,stn5_title + stn5_str]
        }
        
        for key in opts:
            if key != "stations" and  (opts[key] == "y" or opts[key] == "Y"):
                print(strings[key])
            else:
                for i,stn in enumerate(opts[key]):
                    if stn =="Y" or stn == "y":
                        print(strings[key][i])

    def evaluate_turbine(self):
        pass

    def print_states(self):
        print(self.state_1)
        print(self.state_2)
#        print(self.state_3)
        print(self.state_4)
        print(self.state_5)

    def velocity_triangles(self):
        fig, tri = plt.subplots(2)
        tri[0].set_title("Rotor Inlet Velocity Triangles")
        tri[0].arrow(0,0,self.C_theta4,self.C_m4,length_includes_head=True,color="r",width = 0.01,head_width=.75*self.C_theta4/35.2)
        tri[0].arrow(0,self.C_m4,self.C_theta4,0,length_includes_head=True,ls="--",color="r",width = 0.01,head_width=.75*self.C_theta4/35.2)
        tri[0].arrow(0,0,0,self.C_m4,length_includes_head=True,ls="--",color="r",width = 0.01,head_width=.75*self.C_theta4/35.2)

        tri[0].arrow(0,0,self.W_theta4,self.W_m4,length_includes_head=True,color="b",width = 0.01,head_width=.75*self.C_theta4/35.2)
        tri[0].arrow(0,self.W_m4,self.W_theta4,0,length_includes_head=True,ls="--",color="b",width = 0.01,head_width=.75*self.C_theta4/35.2)
        tri[0].arrow(0,0,0,self.W_m4,length_includes_head=True,ls="--",color="b",width = 0.01,head_width=.75*self.C_theta4/35.2)
        c_patch = mpl.lines.Line2D([], [],color='red', label='Absolute Velocity')
        w_patch = mpl.lines.Line2D([], [],color='blue', label='Relative Velocity')
        tri[0].legend(handles=[c_patch, w_patch])
        tri[0].invert_yaxis()
        
        tri[1].set_title("Rotor Outlet Velocity Triangles")
        tri[1].arrow(0,0,0,self.C_m5,length_includes_head=True,color="r",width = 0.01,head_width=0.4)
        tri[1].arrow(0,self.C_m5,0,0,length_includes_head=True,ls="--",color="r",width = 0.01,head_width=0.4)
        tri[1].arrow(0,0,0,self.C_m5,length_includes_head=True,ls="--",color="r",width = 0.01,head_width=0.4)

        tri[1].arrow(0,0,self.W_theta5,self.W_m5,length_includes_head=True,color="b",width = 0.01,head_width=0.4)
        tri[1].arrow(0,self.W_m5,self.W_theta5,0,length_includes_head=True,ls="--",color="b",width = 0.01,head_width=0.4)
        tri[1].arrow(0,0,0,self.W_m5,length_includes_head=True,ls="--",color="b",width = 0.01,head_width=0.4)
        tri[1].invert_yaxis()
        fig.tight_layout()
        for ax in tri.flat:
            ax.set(xlabel='Velocity [m s^-1]', ylabel='Velocity [m s^-1]')
        plt.show()

    def make_hub_and_shroud(self):
        n = 5
        N_p = 20 # Number of points used to make the hub and shroud curves
        """Hub curve"""
        if self.r_4-self.r_h5 > self.z_r:

            R_c = self.z_r
            L = self.r_4 -self.r_h5 - self.z_r
        else:
            R_c = self.r_4-self.r_h5
            L = self.z_r - (self.r_4 - self.r_h5)
        hub_points = np.zeros((N_p,2))
        thetas = np.linspace(3 *np.pi/2, 2 * np.pi, num = N_p)
        for i,theta in enumerate(thetas):
            hub_points[i,0] = R_c * np.cos(theta) 
            hub_points[i,1] = R_c * np.sin(theta) + R_c + self.r_h5

        if L == self.r_4 -self.r_h5 - self.z_r:
            z = hub_points[-1,0]
            r = hub_points[-1,1]
            N_p_L = int(L / R_c * N_p)
            dL = L / N_p_L
            L_mat = np.zeros((N_p_L,2))
            for k in range(N_p_L):
                L_mat[k,0] = z
                L_mat[k,1] = r + (k+1) * dL
            new_hub_points = np.concatenate((hub_points,L_mat),axis=0)
        else: 
            z = hub_points[0,0]
            r = hub_points[0,1]
            N_p_L = int(L / R_c * N_p)
            dL = L / N_p_L 
            L_mat = np.zeros((N_p_L,2))
            for k in range(N_p_L):
                L_mat[N_p_L-1-k,0] =  z - (k+1) * dL
                L_mat[k,1] = r 
            new_hub_points = np.concatenate((L_mat,hub_points),axis=0)
        # print("length added curve:")
        # print(new_hub_points)
        new_new_hub_points = np.zeros((len(new_hub_points),2))
        for q,points in enumerate(new_hub_points):
            new_new_hub_points[q,0] = points[0] + np.abs(new_hub_points[0,0])
            new_new_hub_points[q,1] = points[1]

        hub_points = new_new_hub_points
        plt.plot(hub_points[:,0],hub_points[:,1])

        """Shroud curve"""
        shroud_points = np.zeros((len(hub_points),2))
        z_5 = np.min(hub_points[:,0])
        shroud_z =  np.linspace(z_5, R_c-self.b_4, num = N_p+N_p_L)
        for j,points in enumerate(hub_points):
            zeta = (shroud_z[j] - z_5)/(self.z_r-self.b_4)
            r = self.r_s5 + (self.r_4 - self.r_s5)*zeta**n
            z = shroud_z[j]
            #print(z)
            shroud_points[j,0] = z
            shroud_points[j,1] = r

        """Calculate meridian"""
        m_h = np.zeros((N_p + N_p_L,1))
        m_s = np.zeros((N_p + N_p_L,1))

        for i in range(len(m_h)):
            if i==0:
                m_h[i]=0
                m_s[i]=0
            else:
                m_h[i] = np.sqrt((hub_points[i,0]-hub_points[i-1,0])**2 + (hub_points[i,0]-hub_points[i-1,0])**2) + m_h[i-1]
                m_s[i] = np.sqrt((shroud_points[i,0]-shroud_points[i-1,0])**2 + (shroud_points[i,1]-shroud_points[i-1,1])**2) + m_s[i-1]    
        m_h_mean = np.mean(m_h)    
        m_s_mean = np.mean(m_s)
        print(m_h_mean,m_s_mean)
        plt.plot(shroud_points[:,0],shroud_points[:,1])
        plt.xlim([-.01, 0.06])
        plt.ylim([-.01, 0.06])
        plt.show()

class nozzle:
    def __init__(self,nozzle_inputs,turb):
        self.N_n = turb.n_r + 2 # Number of nozzle blades guess
        self.N_p_n = 1000
        self.turb = turb
        self.sc = nozzle_inputs["sc"]
        self.theta = nozzle_inputs["camber angle"]
        self.ac = nozzle_inputs["ac"]
        self.t_2c = nozzle_inputs["t_2c"]
        self.t_3c = nozzle_inputs["t_3c"]
        self.t_maxc = nozzle_inputs["t_maxc"]
        self.dc = nozzle_inputs["dc"]
        self.radius_ratio = nozzle_inputs["radius ratio"]
        # self.chi_2 = np.arctan((4 * self.bc)/(4*self.ac -1))
        # self.chi_3 = np.arctan()
        self.r_3 = turb.r_3
        self.r_2 = self.r_3 * self.radius_ratio
        self.N_n = turb.n_r + 2 # Number of nozzle blades guess
        self.C_theta3 = turb.r_4 / turb.r_3 * turb.C_theta4
        self.h_03 = turb.state_01.h
        self.P_03 = turb.P_04
        self.b_3 = turb.b_4
        self.state_03 = H2(h=self.h_03,p=self.P_03)
        rho_3_guess = 0
        rho_3_check = self.state_03.rho
        while rho_3_check - rho_3_guess != 0:
            rho_3_guess = rho_3_check
            C_m3 = turb.CFE.mass_flow / (2 * np.pi * self.r_3 * self.b_3 * rho_3_guess)
            C_3 = np.sqrt(C_m3**2 + self.C_theta3**2)
            h_3 = self.h_03 - 0.5 * C_3**2
            s_3 = self.state_03.s
            state_3 = H2(s=s_3,h=h_3)
            rho_3_check = state_3.rho

        self.C_m3 = C_m3
        self.state_3 = state_3
        self.C_3 = C_3
        self.alpha_3_rad = np.arctan(self.C_theta3/self.C_m3)
        self.alpha_3 = self.alpha_3_rad * 180 / np.pi
        self.d_3 = 2 * np.pi * self.r_3 / self.N_n # Pitch is usually s but since s is entropy I am using d
        self.o_3 = self.d_3 * np.sin(self.alpha_3_rad)
        self.c = 1/self.sc * self.d_3

        naca_output = self.calc_naca_profile()
        self.xcs = naca_output[0]
        self.ycs = naca_output[1]
        self.naca_cam = np.array([self.xcs,self.ycs],dtype=object).T
        self.naca_chi = naca_output[2]
        self.naca_suc_surf = naca_output[3]
        self.naca_pres_surf = naca_output[4]
        


    def calc_naca_profile(self):
        N_p_n = self.N_p_n
        check = np.zeros((N_p_n,))
        tcs = np.zeros((N_p_n,))
        suc = np.zeros((N_p_n,2))
        pres = np.zeros((N_p_n,2))

        xcs = np.linspace(0,self.c,N_p_n)
        ycs = np.zeros((N_p_n,))
        # print(np.tan(self.theta))
        # print((self.ac - self.ac**2 - 3/16))
        # print(1 + (4 * np.tan(self.theta))**2 * (self.ac - self.ac**2 - 3/16))
        # print((4 * np.tan(self.theta)))
        bc = (np.sqrt(1 + (4 * np.tan(self.theta))**2 * (self.ac - self.ac**2 - 3/16))-1) / (4 * np.tan(self.theta))
        b = bc*self.c
        a = self.ac * self.c
        # print(b_div_c)
        chi = np.zeros((N_p_n,1))
        chi[0] = np.arctan((4 * b) / (4*a - self.c))
        chi[-1] = np.arctan((4* b) / (3*self.c - 4*a))
        check[0] = chi[0] * 180 / np.pi
        check[-1] = chi[-1] * 180 / np.pi
        for i,xc in enumerate(xcs):
            # print()
            stop = 1000
            k = 0
            yc_guess = 1
            yc_check = 0
            while yc_guess != yc_check:
                k+=1
                yc_guess = yc_check
                yc_num = (xc * (self.c - xc))
                yc_den_1 = (self.c - 2 * a)**2 / (4 * b**2) * yc_guess
                yc_den_2 = (self.c - 2 * a) / (b) * xc
                yc_den_3 = (self.c**2 - 4 * a*self.c) / (4 * b)
                yc_check = yc_num / (yc_den_1 + yc_den_2 - yc_den_3)
                if k == stop:
                    break
            ycs[i] = yc_check
            if i == 1 or i==0:
                pass
            else:
                chi[i-1] = (np.arctan((ycs[i] - ycs[i-2]) / (xcs[i] - xcs[i-2])))
                check[i-1] = chi[i-1] * 180 / np.pi
            if xc <= self.dc*self.c:
                zeta = xc / (self.dc * self.c)
                # print("xc <= dc")

            else:
                zeta = (self.c-xc) / (self.c-self.dc * self.c)
                # print("xc > dc")
            
            # print("zeta:",zeta)
            t_refc = self.t_2c*self.c + (self.t_3c*self.c-self.t_2c*self.c) * xc /(self.c)
            # print("tref:",t_refc)

            # e = np.sqrt(0.4 * self.dc) * (0.03 * (1 - xc/self.c) * (1-zeta) + 0.95)
            e = np.sqrt(0.4 * self.dc) * (0.95 * (1 + xc/self.c) * (1-zeta) +0.05)
            # print("e:",e)

            tcs[i] = t_refc + (self.t_maxc*self.c-t_refc) * zeta **e
            # print("t:",tc)
            if i >= 2:
                pres[i-1,0] = xcs[i-1] + 0.5 * tcs[i-1] * np.sin(chi[i-1])
                pres[i-1,1] = ycs[i-1] - 0.5 * tcs[i-1] * np.cos(chi[i-1])

                suc[i-1,0] = xcs[i-1] - 0.5 * tcs[i-1] * np.sin(chi[i-1])
                suc[i-1,1] = ycs[i-1] + 0.5 * tcs[i-1] * np.cos(chi[i-1])
            elif i==0:
                pres[i,0] = xcs[i] + 0.5 * tcs[i] * np.sin(chi[i])
                pres[i,1] = ycs[i] - 0.5 * tcs[i] * np.cos(chi[i])

                suc[i,0] = xcs[i] - 0.5 * tcs[i] * np.sin(chi[i])
                suc[i,1] = ycs[i] + 0.5 * tcs[i] * np.cos(chi[i])
            
            check[i] = pres[i,1] - suc[i,1]
        if np.sign(chi[-2]) != np.sign(chi[-1]):
            chi[-1] = chi[-1] *-1
        pres[-1,0] = xcs[-1] + 0.5 * tcs[-1] * np.sin(chi[-1])
        pres[-1,1] = ycs[-1] - 0.5 * tcs[-1] * np.cos(chi[-1])

        suc[-1,0] = xcs[-1] - 0.5 * tcs[-1] * np.sin(chi[-1])
        suc[-1,1] = ycs[-1] + 0.5 * tcs[-1] * np.cos(chi[-1])
        print(chi)
        return [xcs,ycs,chi,suc,pres]

    def create_cascade(self):
        rot_cam = np.zeros((self.N_p_n,2))
        rot_suc = np.zeros((self.N_p_n,2))
        rot_pres = np.zeros((self.N_p_n,2))
        blade_2_cam = np.zeros((self.N_p_n,2))
        blade_2_suc = np.zeros((self.N_p_n,2))
        blade_2_pres = np.zeros((self.N_p_n,2))
        gam_guess = 10 / 180 * np.pi
        # r_3c = self.sc * self.N_n / (2 * np.pi)
        c = self.d_3 * 1/self.sc 
        for i,point in enumerate(self.naca_cam):
            suc_point = self.naca_suc_surf[i]
            pres_point = self.naca_pres_surf[i]
            rotated_camber = rotate_blade_point(point,gam_guess,c,self.r_3)
            rotated_suc = rotate_blade_point(suc_point,gam_guess,c,self.r_3)
            rotated_pres = rotate_blade_point(pres_point,gam_guess,c,self.r_3)

            rot_cam[i] = rotated_camber[0]
            rot_suc[i] = rotated_suc[0]
            rot_pres[i] = rotated_pres[0]

            blade_2_camber = make_second_blade(rotated_camber,self.N_n,i)
            blade_2_suction = make_second_blade(rotated_suc,self.N_n,i)
            blade_2_pressure = make_second_blade(rotated_pres,self.N_n,i)

            blade_2_cam[i,0] = blade_2_camber[0]
            if i > 0 and blade_2_cam[i-1,0] > blade_2_cam[i,0]:
                blade_2_cam[i-1,0] = - blade_2_cam[i-1,0]
            blade_2_cam[i,1] = blade_2_camber[1]

            blade_2_pres[i,0] = blade_2_pressure[0]
            if blade_2_cam[i,0] < np.abs*(blade_2_pres[i,0]):
                blade_2_pres[i,0] = - blade_2_pres[i,0]
            blade_2_pres[i,1] = blade_2_pressure[1]

            blade_2_suc[i,0] = blade_2_suction[0]
            if blade_2_cam[i,0] > np.abs(blade_2_suc[i,0]):
                blade_2_suc[i,0] = - blade_2_suc[i,0]
            blade_2_suc[i,1] = blade_2_suction[1]

        thetas = np.linspace(0, 2 * np.pi, num = 100)
        circ = np.zeros((100,2))
        for i,theta in enumerate(thetas):
            circ[i,0] = self.r_3 * (np.cos(theta))
            circ[i,1] = self.r_3 * (np.sin(theta))
        plt.plot(circ[:,0],circ[:,1])
        plt.plot(blade_2_suc[:,0],blade_2_suc[:,1])
        plt.plot(blade_2_pres[:,0],blade_2_pres[:,1],marker=".")
        plt.plot(blade_2_cam[:,0],blade_2_cam[:,1])
        plt.plot(rot_cam[:,0],rot_cam[:,1])
        plt.plot(rot_suc[:,0],rot_suc[:,1])
        plt.plot(rot_pres[:,0],rot_pres[:,1])
        plt.show()

    def find_setting_angle(self):
        pass

    def plot_naca_norm(self):
        # FIXME - Make this defined within the nozzle class
        # nozx = noz_out[0]
        # nozy = noz_out[1]
        # nozchi = noz_out[2]
        # nozsuc = noz_out[3]
        # nozpres = noz_out[4]
        # print(nozpres)
        # print(nozsuc)
        # r_2 = noz.t_2c/2
        # r_3 = noz.t_3c/2
        # thetas = np.linspace(0, 2 * np.pi, num = 100)
        # le = np.zeros((100,2))
        # te = np.zeros((100,2))
        # for i,theta in enumerate(thetas):
        #     le[i,0] = r_2 * (np.cos(theta))
        #     le[i,1] = r_2 * (np.sin(theta))
        #     te[i,0] = r_3 * (np.cos(theta)) + nozx[-1]
        #     te[i,1] = r_3 * (np.sin(theta)) + nozy[-1]
        # plt.plot(le[:,0],le[:,1])
        # plt.plot(te[:,0],te[:,1])
        # plt.plot(nozx,nozy)
        # plt.plot(nozpres[:,0],nozpres[:,1],marker = ".")
        # plt.plot(nozsuc[:,0],nozsuc[:,1],marker = ".")
        # plt.xlim([-0.01,1.01])
        # plt.ylim([-0.1,0.1])
        # plt.show()
        pass
def rotate_blade_point(point,gam_guess,c,r_3):
    rot_point = np.zeros((2,))
    rot_point[0] = (point[0] - c) * np.cos(gam_guess) + point[1] * np.sin(gam_guess)
    rot_point[1] = r_3 - (point[0] - c) * np.sin(gam_guess) + point[1] * np.cos(gam_guess)
    theta1 = np.arctan(rot_point[0]/rot_point[1])
    rb1 = np.sqrt(rot_point[0]**2 + rot_point[1]**2)

    return [rot_point,theta1,rb1]

def make_second_blade(rot_point,N_n,i):
    blade_2 = np.zeros((2,))
    rb2 = rot_point[2]
    theta2 = rot_point[1] + 2 * np.pi / N_n
    tan = np.tan(theta2)

    blade_2[0] = np.sqrt(((rb2*tan)**2) / (1 + tan**2))

    blade_2[1] = np.sqrt(rb2**2 - blade_2[0]**2)
    return blade_2

def find_turb(init_cfe,init_turb):
    new_turb = init_turb
    new_cfe = init_cfe
    i = new_turb.i
    prt_opt ="N"
    while np.round(new_turb.eta_ts,4) != np.round(new_turb.eta_ts_loss,4):
        turb = new_turb
        cfe = new_cfe
        if prt_opt =="Y":
            print("Turbine iteration:",turb.i)
            print("Turbine efficiency:",turb.eta_ts)
            print("Turbine specific speed:",turb.N_s)
            print("Turbine specific diameter:",turb.D_s)
            print("N_s * D_s:",turb.N_s*turb.D_s)
            print("Turbine total-to-total enthalpy drop:",turb.deltah_0)
            print("Turbine total-to-static pressure ratio:",turb.PR_ts)
            print("Turbine enthalpy losses:",turb.h_loss)
        i += 1
        new_h_5ss = turb.h_01 - turb.deltah_0 / turb.eta_ts
        turb.state_5ss = H2(h=new_h_5ss,s=turb.state_01.s)
        new_P_5 =turb.state_5ss.p
        new_PR = turb.state_01.p / new_P_5

        new_dynamic_turb_inputs = {
            "PR_ts" : new_PR,
            "eta_ts" : turb.eta_ts_loss,
            "h_0ss" : 0,
        }

        new_cfe = CFE(cfe.static_cfe_inputs,new_dynamic_turb_inputs,i)
        new_delta_h0 = new_cfe.work_rate / new_cfe.mass_flow
        new_delta_h0ss = new_delta_h0 / turb.eta_ts_loss

        new_dynamic_turb_inputs = {
            "PR_ts" : new_PR,
            "eta_ts" : turb.eta_ts_loss,
            "h_0ss" : new_delta_h0ss
            }
        new_static_turb_inputs = turb.static_turb_inputs
        new_turb = turbine(new_cfe,new_static_turb_inputs,new_dynamic_turb_inputs,i)
    return new_turb
            
if __name__ == "__main__":
    pass
    # working_fluid_inputs = {
    #     "mu" : "-0.00000000000144 *T**2 + 0.0000000169 *T+ 0.00000464",#[Pa-s] - Kinematic Viscosity
    # }

    # nozzle_inputs = {
    #     "radius ratio" : 1.1,
    #     "camber angle" : 0,
    #     "ac" : 0.7,
    #     "t_2c" : 0.025,
    #     "t_3c" : 0.012,
    #     "t_maxc" : 0.06,
    #     "dc" : 0.4,
    #     "sc" : 0.75
    # }
    # static_cfe_inputs = {
    #     "inner_radius" : 0.056, #Channel inner radius [m]
    #     "outer_radius" : 0.066, #Channel outer radius [m]
    #     "length" : 0.94, #CFE channel length [m]
    #     "rpm" : 7000, #CFE inner SiC cylinder revolutions per minute
    #     "mass_flow" : 0.108, #CFE channel mass flow rate [kg s^-1]

    #     "temp" : 450, #[K]
    #     "press" : 12.5 #MPa - Turbine Inlet Pressure
    # } 

    # dynamic_turb_inputs = {
    #     "PR_ts" : 4,
    #     "eta_ts" : 0.9,
    #     "h_0ss" : 0,
    #     "N_s" : 0
    # }
    # opts = {
    #     "dim" : "Y",
    #     "prelim" : "y",
    #     "geom" : "y",
    #     "stations" : ["Y","y"]
    # }
    # test_cfe = CFE(static_cfe_inputs,dynamic_turb_inputs,1)

    # init_turb = turbine(test_cfe,test_cfe.static_turb_inputs,dynamic_turb_inputs,1)

    # test_turb = find_turb(test_cfe,init_turb)
    # # test_turb.print_turbine(opts)
    # # test_turb.make_hub_and_shroud()
    # # test_turb.velocity_triangles()
    # # test_turb.print_states()
    # # """"prelim calcs"""

    # """Checks:"""
    # state_5ss = H2(s = test_turb.state_01.s,p = test_turb.state_5.p)
    # eff_check = (test_turb.h_01 - test_turb.h_05) / (test_turb.h_01 - test_turb.state_5ss.h)
    # print(test_turb.h_01 - test_turb.h_05)
    # print(test_turb.deltah_0)
    # print(eff_check)
    # print()
    # print("Mach No. Check:")
    # M_inlet = test_turb.C_4/test_turb.state_4.a
    # M_outlet = test_turb.C_5/test_turb.state_5.a
    # print("Inlet Mach No.:",M_inlet)
    # print("Outlet Mach No.",M_outlet)

    # print()
    # print("Inlet Relative Flow Angle:")
    # W_theta4 = test_turb.C_theta4-test_turb.U_4
    # W_m4=test_turb.C_m4
    # beta_4 = np.arctan(W_theta4/W_m4)
    # print(beta_4*180/np.pi)

    # print()
    # print("Axial Length Check:")
    # print(1.5*test_turb.b_4,"<?=",test_turb.z_r)

    # print()
    # print("Outlet Meridional Velocity and Shroud Ratios:")
    # vel_ratio1 = test_turb.C_m5/test_turb.U_4
    # rad_ratio = test_turb.r_s5/test_turb.r_4
    # print(vel_ratio1)
    # print(rad_ratio)
    

    # print()
    # print("Meridional Velocity Ratio:")
    # MVR = test_turb.C_m5/test_turb.C_m4
    # print(MVR)
    
    # print()
    # print("Stage Reaction:")
    # R = (test_turb.state_4.h-test_turb.state_5.h)/(test_turb.state_01.h-test_turb.state_05.h)
    # print(R)
