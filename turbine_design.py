import numpy as np
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
    def __init__(self,static_cfe_inputs,dynamic_turb_inputs,H2):

        """ Entry channel fluid properties """
        self.T_in = static_cfe_inputs["temp"] #Temperature [K]
        self.P_in = static_cfe_inputs["press"]*dynamic_turb_inputs["PR_ts_guess"] * 1e6
        self.P_out = static_cfe_inputs["press"] * 1e6
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
        self.work_rate = self.calc_work_rate() + 55#Power requirement for turbine [W]
        print(self.work_rate)
        st_turb_inputs = self.calc_static_turb_inputs()
        self.static_turb_inputs = {
            "T_01" : st_turb_inputs[0],
            "C_m1" : st_turb_inputs[1],
            "C_theta1" : st_turb_inputs[2]
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
            G = 10**(-0.006360 * (np.log10(Re))**3 + 0.1349 * (np.log10(Re))**2 + 0.8850 * (np.log10(Re)) +1.610 )

        eta =15.999/22.085 
    
        G_lam = 4*np.pi * eta /( (1-eta)*(1-eta**2))*Re
        Nu_omega = G/G_lam
        G_lam = 2 * self.eta /( (1-self.eta)*(1-self.eta**2))*Re
        G = Nu_omega * G_lam
        T_1 = G * 2*np.pi*self.L*self.rho*self.nu**2
        T_lam = 4*np.pi*self.L*self.mu*self.omega/(self.R_i**-2 - self.R_o**-2)
        T_2 = T_lam * Nu_omega 
        work = T_1*self.omega
        return work

    def calc_static_turb_inputs(self):
        U_thetaB = self.omega * self.R_i**2 / (self.R_o + self.R_i)
        U_mB = self.mass_flow / self.cfe_state.rho/self.annulus_area
        U_B = np.sqrt(U_thetaB**2 + U_mB**2)
        M_B = U_B/np.sqrt(self.cfe_state.gamma * 4124.2 * self.cfe_state.t)
        T_01 = self.T_in * ( 1 + (self.cfe_state.gamma - 1) / 2 * M_B**2 )
        return [T_01,U_thetaB,U_mB]

def eval_scalar_input(inpt):
    if isinstance(inpt, str):
        return eval(inpt)
    else:
        return inpt 

class station:
    def __init__(self,turb):
        pass
class whitfield_turbine:
    def __init__(self,CFE,static_turb_inputs,dynamic_turb_inputs,i,H2):
        self.CFE = CFE
        self.i = i
        self.state_01 = H2(t=static_turb_inputs["T_01"], p = self.CFE.P_in)
        self.T_01 = static_turb_inputs["T_01"]
        self.h_01 = self.state_01.h
        self.delta_W = self.CFE.work_rate / self.CFE.mass_flow
        self.S = self.delta_W / self.h_01 #Power ratio
        self.T_03 = (1 - self.S) * self.T_01


class turbine:
    def __init__(self,CFE,static_turb_inputs,dynamic_turb_inputs,i,H2):
        self.CFE = CFE
        self.i = i
        #self.number_stations = static_turb_inputs["number of stations"]
        self.stations = self.make_stations()
        
        #self.eta_tt_guess = dynamic_turb_inputs["total to total efficiency guess"]
        self.eta_ts_guess = dynamic_turb_inputs["eta_ts"]
        self.PR_ts_guess = dynamic_turb_inputs["PR_ts_guess"]

        self.deltah_0 = self.CFE.work_rate/self.CFE.mass_flow
        self.state_01 = H2(t=static_turb_inputs["T_01"], p = self.CFE.P_in)
        self.h_05 = self.state_01.h - self.deltah_0
        self.state_05 = H2(h = self.h_05, P=self.CFE.P_out)
        self.epsilon_b = 0.4e-3
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
            self.N_s = dynamic_turb_inputs["N_s"]

        """Preliminary Calcs: Rotor"""
        """Isentropic Values"""
        self.v_s = 0.667 #0.737 * self.N_s**0.2
        self.C_0 = np.sqrt(2 * self.h_0ss)
        """Total and static states at station 4"""
        self.U_4 = self.C_0 * self.v_s
        self.r_4 = self.U_4/self.CFE.omega
        self.alpha_4 = 90 - 1* (10.8 + 14.2*self.N_s**2)
        self.alpha_4_rad = self.alpha_4 / 180 * np.pi
        self.P_04 = self.CFE.P_in - self.state_01.rho * self.h_0ss*(1-self.eta_ts)/4
        self.C_theta4 = self.U_4 * self.eta_ts / (2*self.v_s**2) #??? is this okay????
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
        self.a_4 = np.sqrt(self.state_4.gamma * 4124.2 * self.state_4.T)
        self.M_4 = self.C_4 / self.a_4

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
        self.a_5 = np.sqrt(self.state_5.gamma * 4124.2 * self.state_5.T)

        

        self.z_r = 1.5 * (self.b_5)
        self.n_r = np.round(np.pi/30*(110-self.alpha_4)*np.tan(self.alpha_4*np.pi/180)) #Glassman 1972
        self.q_5 = 2 * np.pi * self.r_5 / self.n_r
        self.o_5 = self.q_5 * self.C_m5 / self.W_5
        one_over_PR = (1 - (self.C_0**2 / (2*self.state_01.cp*self.state_01.t)))**(self.state_01.gamma/(self.state_01.gamma-1))
        self.PR = 1/one_over_PR
        self.r_3 = self.r_4 + 2 * self.b_4 * np.cos(self.alpha_4*np.pi/180)
        self.t_lead = 0.04 * self.r_4
        self.t_trail = 0.02 * self.r_4


        """Preliminary Calcs: The Bowl"""
        self.state_1 = H2(t=self.CFE.T_in,p=self.CFE.P_in)
        self.r_1 = self.CFE.R_i + self.CFE.h / 2
        self.A_1 = self.CFE.annulus_area
        self.C_m1 = static_turb_inputs["C_m1"]
        self.C_theta1 = static_turb_inputs["C_theta1"]
        self.C_1 = np.sqrt(self.C_m1**2 + self.C_theta1**2)
        self.h_01 = self.state_1.h + 0.5 * self.C_1**2
        
        self.h_02 = self.h_01
        self.A_2 = 2 * np.pi * self.r_4 * self.b_4
        self.C_m2 = self.C_m1 * self.A_1 / self.A_2
        print("C_m2:",self.C_m2)
        self.r_2 = 1.1 * self.r_3
        print(self.r_3)
        print(self.r_2)
        self.C_theta2 = self.C_theta1 * self.r_1/self.r_2
        print(self.C_theta2)
        self.alpha_2 = np.tan(self.C_theta2/self.C_m2)*180/np.pi
        print(self.alpha_2)
        self.C_2 = np.sqrt(self.C_m2**2 + self.C_theta2**2)
        self.h_2 = self.h_02 - 0.5 * self.C_2**2
        self.s_2 = self.state_1.s
        self.state_2 = H2(h=self.h_2,s=self.s_2)
        C_m2_check = self.CFE.mass_flow/(2 * np.pi * self.r_2 * self.state_2.rho * self.b_4)
        print("Check:",C_m2_check)


        """Preliminary Calcs: Nozzle"""


    def nozzle_design(self,s_div_c = 0.75, theta = 0, a_div_c = 0.5, t_2_div_c = 0.025,
        t_3_div_c = 0.012, t_max_div_c = 0.06, d_div_c = 0.4,geom_gamma_3 = 20):
        pass
        

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
        print("Calculated eta:", eta_ts)
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

        M_5_rel = self.W_5 /self.a_5

        """Incidence Losses"""
        L_i_a = (0.75 / (2 * (1 - self.M_4**2)))
        L_i_b = ((self.C_m4 / self.C_m5) * (self.C_m4 / self.U_4) * np.tan(self.alpha_4*np.pi/180) - 1)**2
        L_i =  L_i_a * L_i_b * self.U_4**2
        print("Incidence Losses:",L_i)

        """Passage Losses"""
        L_p_f = L_h/D_h
        L_p_sf = 0.68 * (1 - (r_t/self.r_4)**2) * (np.cos(beta_t) / (b_t/chord_len))
        L_p = m_f * 0.11 * (L_p_f + L_p_sf) * (self.W_4**2 + W_t**2)/2
        print("Passage Losses:",L_p)

        """Windage Losses"""
        L_w = K_f * (rho_avg * self.U_4**3 * self.r_4**2) / (2 * self.CFE.mass_flow * self.W_5**2)
        print("Windage Losses:",L_w)

        """Trailing Edge Losses"""
        g = 9.81
        dP_part2 = (self.n_r * self.t_lead / (np.pi * (self.r_h5 + self.r_s5) * np.cos(self.beta_5_rad)))**2
        deltaP_0_rel = self.state_5.rho * self.W_5**2 / (2 * g) * dP_part2
        
        L_t_den = self.state_5.p * (1 + (self.W_5**2/(2 * self.state_5.t * self.state_5.cp)))**(self.state_5.gamma/(self.state_5.gamma-1))
        L_t = 2 /(self.state_5.gamma * M_5_rel**2 )  * (deltaP_0_rel/L_t_den)
        print("Trailing edge losses:",L_t)

        """Exit Energy Loss"""
        L_e = self.C_5**2 / 2
        print("Exit Energy Losses:",L_e)

        """Nozzle Losses"""

        print("Nozzle Losses:",L_n)

        h_losses = ( L_e + L_i  + L_p + L_w + L_t + L_n)
        return h_losses


    def print_turbine(self):
        turb_title = f'Turbine iteration: {self.i}\n\n'
        stn1_str_title = f'Station 1:\n\n'
        stn2_str_title = f'Station 2:\n\n'
        stn3_str_title = f'Station 3:\n\n'
        stn4_str_title = f'Station 4:\n\n'
        stn5_str_title = f'Station 5:\n\n'
        stn6_str_title = f'Station 6:\n\n'
        isen_vals_str = f' \
Preliminary calculated total to static efficicency: {self.eta_ts}\n \
Isentropic enthalpy drop: {self.h_0ss} [kJ kg^-1]\n \
Total-total enthalpy drop: {self.deltah_0} [kJ kg^-1]\n \
Total-to-static pressure ratio: {self.PR}\n \
Specific speed: {self.N_s}\n \
Turbine velocity ratio: {self.v_s}\n \
Spouting Velocity: {self.C_0} [m s^-1]\n \
Turbine rotor inlet blade speed: {self.U_4} [m s^-1]\n\n'

        geom_str = f' \
Rotor inlet radius: {self.r_4*1000} [mm]\n \
Rotor inlet blade height: {self.b_4*1000} [mm]\n \
Rotor outlet hub radius: {self.r_h5 * 1000} [mm]\n \
Rotor outlet shroud radius: {self.r_s5 * 1000} [mm]\n \
Rotor outlet blade height: {self.b_5 * 1000} [mm]\n \
Rotor axial length: {self.z_r * 1000} [mm]\n \
Number of rotor blades: {self.n_r}\n \
Rotor inlet relative flow angle: {self.beta_4} [degrees]\n\n'

        stn4_str = f' \
Rotor inlet stagnation pressure: {self.P_04/1000} [kPa]\n \
Rotor inlet radius: {self.r_4*1000} [mm]\n \
Absolute meridional velocity: {self.C_m4} [m s^-1]\n \
Absolute tangential velocity: {self.C_theta4} [m s^-1]\n \
Absolute velocity: {self.C_4} [m s^-1]\n \
Rotor inlet absolute flow angle: {self.alpha_4} [degrees]\n \
Rotor inlet relative flow angle: {self.beta_4} [degrees]\n\n'
        print(turb_title)
        print(isen_vals_str + geom_str + stn4_str_title+stn4_str)

    def calculate_losses(self):
        pass

    def evaluate_turbine(self):
        pass

    def print_states(self):
        print(self.state_1)
        print(self.state_2)
#        print(self.state_3)
        print(self.state_4)
        print(self.state_5)
def find_turb(init_cfe,init_turb):
    pass

            
if __name__ == "__main__":
    H2 = Fluid("Hydrogen", prop_files)

    working_fluid_inputs = {
        "mu" : "-0.00000000000144 *T**2 + 0.0000000169 *T+ 0.00000464",#[Pa-s] - Kinematic Viscosity
    }

    static_cfe_inputs = {
        "inner_radius" : 0.056, #Channel inner radius [m]
        "outer_radius" : 0.058, #Channel outer radius [m]
        "length" : 0.84, #CFE channel length [m]
        "rpm" : 7000, #CFE inner SiC cylinder revolutions per minute
        "mass_flow" : 0.108, #CFE channel mass flow rate [kg s^-1]

        "temp" : 450, #[K]
        "press" : 12.5 #MPa - Turbine Outlet Pressure
    } 


    dynamic_turb_inputs = {
        "PR_ts_guess" : 1.00096,
        "eta_ts" : 0.9,
        "h_0ss" : 0,
        "N_s" : 0
    }


    test_cfe = CFE(static_cfe_inputs,dynamic_turb_inputs,H2)

    test_turb = turbine(test_cfe,test_cfe.static_turb_inputs,dynamic_turb_inputs,1,H2)
    test_turb.print_turbine()
    test_turb.calc_losses()
    test_turb.print_states()
    # """"prelim calcs"""

    """Checks:"""
    print()
    print("Mach No. Check:")
    M_inlet = test_turb.C_4/np.sqrt(test_turb.state_4.gamma*4124.2*test_turb.state_4.t)
    M_outlet = test_turb.C_5/np.sqrt(test_turb.state_5.gamma*4124.2*test_turb.state_5.t)
    print(test_turb.state_4.a)
    print(M_inlet,M_outlet)

    print()
    print("Inlet Relative Flow Angle:")
    W_theta4 = test_turb.C_theta4-test_turb.U_4
    W_m4=test_turb.C_m4
    beta_4 = np.arctan(W_theta4/W_m4)
    print(beta_4*180/np.pi)

    print()
    print("Axial Length Check:")
    print(1.5*test_turb.b_4,"<?=",test_turb.z_r)

    print()
    print("Outlet Meridional Velocity and Shroud Ratios:")
    vel_ratio1 = test_turb.C_m5/test_turb.U_4
    rad_ratio = test_turb.r_s5/test_turb.r_4
    print(vel_ratio1)
    print(rad_ratio)

    print()
    print("Meridional Velocity Ratio:")
    MVR = test_turb.C_m5/test_turb.C_m4
    print(MVR)
    
    print()
    print("Stage Reaction:")
    R = (test_turb.state_4.h-test_turb.state_5.h)/(test_turb.state_01.h-test_turb.state_05.h)
    print(R)

    # PR = (1 - (C_0**2 / (2*test_turb.state_01.cp*test_turb.state_01.t)))**(test_turb.state_01.gamma/(test_turb.state_01.gamma-1))
    # print(1/PR)
    # print(state_5.gamma)
    # print(test_turb.state_01.gamma)
    # print(0.02*r_4)

    # q_5 = 2 * np.pi * r_5 / n_r
    # print("Rotor outlet mean blade pitch:",q_5)

    # W_m5 = C_m5
    # W_theta5 = -(r_5 * test_turb.CFE.omega)
    # W_5 = np.sqrt(W_m5**2 + W_theta5**2)
    # o_t = q_5 * C_m5 /W_5
    # print("Mean throat width:",o_t)


    # print("Nozzle Design")
    # b_2 = b_4
    # b_3 = b_2
    # deltaR = 2 * b_4 * np.cos(alpha_4*np.pi/180)
    # print(deltaR)

    # r_3 = r_4 + deltaR
    # print("Nozzle outlet radius:",r_3)

    # C_theta3 = C_theta4 * (r_4/r_3)
    # print("Nozzle outlet tangential velocity:",C_theta3)

    # C_m3 = C_m4 / np.tan(alpha_4)
    # print("Nozzle outlet meridional velocity:",C_m3)
    









# # RPMs
#     rpms = [1000]
#     for i in range(500):
#         rpms.append(rpms[i]+20)
    
#     rpms_load = []
#     rpms_flow = []
#     rpms_N_s = []
#     rpms_D_s = []
#     for rpm in rpms:
#         Turb = turb(T,P,R_i,R_o,L,mdot,rpm)
#         rpms_load.append(Turb.load_coef)
#         rpms_flow.append(Turb.flow_coef)
#         rpms_N_s.append(Turb.N_s)
#         rpms_D_s.append(Turb.D_s)
        
#     plt.title("Loading Coefficient vs RPMs")
#     plt.xlabel("RPMs (s^-1)")
#     plt.ylabel("Psi")
#     plt.plot(rpms,rpms_load)
#     plt.show()
    
#     plt.title("Flow Coefficient vs RPMs")
#     plt.xlabel("RPMs (s^-1)")
#     plt.ylabel("Phi")
#     plt.plot(rpms,rpms_flow)
#     plt.show()
    
#     plt.title("Specific Speed vs RPMs")
#     plt.xlabel("RPMs (s^-1)")
#     plt.ylabel("N_s")
#     plt.plot(rpms,rpms_N_s)
#     plt.show()
    
#     plt.title("Specific Diameter vs RPMs")
#     plt.xlabel("RPMs (s^-1)")
#     plt.ylabel("D_s")
#     plt.plot(rpms,rpms_D_s)
#     plt.show()
        
# #Increasing Outer Radius
#     T = 300
#     P = 11e6
#     R_i = 0.049
#     R_o = 0.051
#     rpm = 7000
#     L = 0.84
#     mdot = 0.108
    
#     R_os = [0.051]
#     for i in range(500):
#         R_os.append(R_os[i]+.005)
    
#     R_os_load = []
#     R_os_flow = []
#     R_os_N_s = []
#     R_os_D_s = []
#     for R_o in R_os:
#         Turb = turb(T,P,R_i,R_o,L,mdot,rpm)
#         R_os_load.append(Turb.load_coef)
#         R_os_flow.append(Turb.flow_coef)
#         R_os_N_s.append(Turb.N_s)
#         R_os_D_s.append(Turb.D_s)
        
#     plt.title("Loading Coefficient vs Outer Radius")
#     plt.xlabel("Outer Radius (m)")
#     plt.ylabel("Psi")
#     plt.plot(R_os,R_os_load)
#     plt.show()
    
#     plt.title("Flow Coefficient vs Outer Radius")
#     plt.xlabel("Outer Radius (m)")
#     plt.ylabel("Phi")
#     plt.plot(R_os,R_os_flow)
#     plt.show()
    
#     plt.title("Specific Speed vs Outer Radius")
#     plt.xlabel("Outer Radius (m)")
#     plt.ylabel("N_s")
#     plt.plot(R_os,R_os_N_s)
#     plt.show()
    
#     plt.title("Specific Diameter vs Outer Radius")
#     plt.xlabel("Outer Radius (m)")
#     plt.ylabel("D_s")
#     plt.plot(R_os,R_os_D_s)
#     plt.show()
    
#     #Decreasing outer radius w/o changing channel height
#     T = 300
#     P = 11e6
#     R_i = 0.049
#     R_o = 0.051
#     rpm = 7000
#     L = 0.84
#     mdot = 0.108
    
#     R_os = [0.051]
#     R_is = [0.049]
#     for i in range(90):
#         R_os.append(R_os[i]-0.00005)
#         R_is.append(R_is[i]-0.0005)
    
#     R_os_load = []
#     R_os_flow = []
#     R_os_N_s = []
#     R_os_D_s = []
#     for i,R_o in enumerate(R_os):
#         Turb = turb(T,P,R_is[i],R_o,L,mdot,rpm)
#         R_os_load.append(Turb.load_coef)
#         R_os_flow.append(Turb.flow_coef)
#         R_os_N_s.append(Turb.N_s)
#         R_os_D_s.append(Turb.D_s)
        
#     plt.title("Loading Coefficient vs Outer Radius")
#     plt.xlabel("Outer Radius (m)")
#     plt.ylabel("Psi")
#     plt.plot(R_os,R_os_load)
#     plt.show()
    
#     plt.title("Flow Coefficient vs Outer Radius")
#     plt.xlabel("Outer Radius (m)")
#     plt.ylabel("Phi")
#     plt.plot(R_os,R_os_flow)
#     plt.show()
    
#     plt.title("Specific Speed vs Outer Radius")
#     plt.xlabel("Outer Radius (m)")
#     plt.ylabel("N_s")
#     plt.plot(R_os,R_os_N_s)
#     plt.show()
    
#     plt.title("Specific Diameter vs Outer Radius")
#     plt.xlabel("Outer Radius (m)")
#     plt.ylabel("D_s")
#     plt.plot(R_os,R_os_D_s)
#     plt.show()
    
#     #Increasing mass flow rate
#     T = 300
#     P = 11e6
#     R_i = 0.049
#     R_o = 0.051
#     rpm = 7000
#     L = 0.84
#     mdot = 0.108
    
#     mdots = [0.108]
#     for i in range(90):
#         mdots.append(mdots[i]+.05)
    
#     mdots_load = []
#     mdots_flow = []
#     mdots_N_s = []
#     mdots_D_s = []
#     for i,mdot in enumerate(mdots):
#         Turb = turb(T,P,R_i,R_o,L,mdot,rpm)
#         mdots_load.append(Turb.load_coef)
#         mdots_flow.append(Turb.flow_coef)
#         mdots_N_s.append(Turb.N_s)
#         mdots_D_s.append(Turb.D_s)
        
#     plt.title("Loading Coefficient vs Mass Flow Rate")
#     plt.xlabel("Mass Flow (kg/s)")
#     plt.ylabel("{\Psi}")
#     plt.plot(mdots,mdots_load)
#     plt.show()
    
#     plt.title("Flow Coefficient vs Mass Flow Rate")
#     plt.xlabel("Mass Flow (kg/s)")
#     plt.ylabel("{\Phi}")
#     plt.plot(mdots,mdots_flow)
#     plt.show()
    
#     plt.title("Specific Speed vs Mass Flow Rate")
#     plt.xlabel("Mass Flow (kg/s)")
#     plt.ylabel("N_s")
#     plt.plot(mdots,mdots_N_s)
#     plt.show()
    
#     plt.title("Specific Diameter vs Mass Flow Rate")
#     plt.xlabel("Mass Flow (kg/s)")
#     plt.ylabel("D_s")
#     plt.plot(mdots,mdots_D_s)
#     plt.show()
    
#         #Increasing Temp
#     T = 300
#     P = 11e6
#     R_i = 0.049
#     R_o = 0.051
#     rpm = 7000
#     L = 0.84
#     mdot = 0.108
    
#     mdots = [273]
#     for i in range(527):
#         mdots.append(mdots[i]+1)
    
#     mdots_load = []
#     mdots_flow = []
#     mdots_N_s = []
#     mdots_D_s = []
#     for i,T in enumerate(mdots):
#         Turb = turb(T,P,R_i,R_o,L,mdot,rpm)
#         mdots_load.append(Turb.load_coef)
#         mdots_flow.append(Turb.flow_coef)
#         mdots_N_s.append(Turb.N_s)
#         mdots_D_s.append(Turb.D_s)
        
#     plt.title("Loading Coefficient vs Temperature")
#     plt.xlabel("Temperature (K)")
#     plt.ylabel("Psi")
#     plt.plot(mdots,mdots_load)
#     plt.show()
    
#     plt.title("Flow Coefficient vs Temperature")
#     plt.xlabel("Temperature (K)")
#     plt.ylabel("Phi")
#     plt.plot(mdots,mdots_flow)
#     plt.show()
    
#     plt.title("Specific Speed vs Temperature")
#     plt.xlabel("Temperature (K)")
#     plt.ylabel("N_s")
#     plt.plot(mdots,mdots_N_s)
#     plt.show()
    
#     plt.title("Specific Diameter vs Temperature")
#     plt.xlabel("Temperature (K)")
#     plt.ylabel("D_s")
#     plt.plot(mdots,mdots_D_s)
#     plt.show()
    
#             #Increasing Pressure
#     T = 300
#     P = 11e6
#     R_i = 0.049
#     R_o = 0.051
#     rpm = 7000
#     L = 0.84
#     mdot = 0.108
    
#     mdots = [0.1]
#     for i in range(527):
#         mdots.append(mdots[i]+.05)
    
#     mdots_load = []
#     mdots_flow = []
#     mdots_N_s = []
#     mdots_D_s = []
#     for i,P in enumerate(mdots):
#         Turb = turb(T,P*1e6,R_i,R_o,L,mdot,rpm)
#         mdots_load.append(Turb.load_coef)
#         mdots_flow.append(Turb.flow_coef)
#         mdots_N_s.append(Turb.N_s)
#         mdots_D_s.append(Turb.D_s)
        
#     plt.title("Loading Coefficient vs Pressure")
#     plt.xlabel("Pressure (MPa)")
#     plt.ylabel("Psi")
#     plt.plot(mdots,mdots_load)
#     plt.show()
    
#     plt.title("Flow Coefficient vs Pressure")
#     plt.xlabel("Pressure (MPa)")
#     plt.ylabel("Phi")
#     plt.plot(mdots,mdots_flow)
#     plt.show()
    
#     plt.title("Specific Speed vs Pressure")
#     plt.xlabel("Pressuer (MPa)")
#     plt.ylabel("N_s")
#     plt.plot(mdots,mdots_N_s)
#     plt.show()
    
#     plt.title("Specific Diameter vs Pressure")
#     plt.xlabel("Pressure (MPa)")
#     plt.ylabel("D_s")
#     plt.plot(mdots,mdots_D_s)
#     plt.show()