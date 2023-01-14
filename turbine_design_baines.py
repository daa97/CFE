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
        """ Entry channel fluid properties """
        self.T_in = static_cfe_inputs["temp"] #Temperature [K]
        self.P_in = static_cfe_inputs["press"]*dynamic_turb_inputs["PR_ts"] * 1e6
        print("CFE Inlet Pressure:",self.P_in/1e6,"[MPa]")
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
        self.work_rate = self.calc_work_rate()  #Power requirement for turbine [W]
        print("Work Rate:",self.work_rate,"[W]")
        st_turb_inputs = self.calc_static_turb_inputs()
        # print(55 *static_cfe_inputs["rpm"]/7000) # checking the bearing resistance relative to rpm
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
        T_bearings = 0.072647
        # print("Viscous torque:",T_1)
        T_1 += T_bearings
        work = T_1*self.omega
        return work

    def calc_static_turb_inputs(self):
        U_thetaB = self.omega * self.R_i**2 / (self.R_o + self.R_i)
        U_mB = self.mass_flow / self.cfe_state.rho/self.annulus_area
        U_B = np.sqrt(U_thetaB**2 + U_mB**2)
        M_B = U_B/self.cfe_state.a
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
        #self.number_stations = static_turb_inputs["number of stations"]
        self.stations = self.make_stations()
        
        #self.eta_tt_guess = dynamic_turb_inputs["total to total efficiency guess"]
        self.eta_ts_guess = dynamic_turb_inputs["eta_ts"]
        self.PR_ts = dynamic_turb_inputs["PR_ts"]

        self.deltah_0 = self.CFE.work_rate/self.CFE.mass_flow
        self.T_01 = static_turb_inputs["T_01"]
        self.state_01 = H2(t=self.T_01, p = self.CFE.P_in)
        self.h_05 = self.state_01.h - self.deltah_0
        self.state_05 = H2(h = self.h_05, P=self.CFE.P_out)
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
        self.v_s = 0.696 # 0.737 * self.N_s**0.2
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
        self.a_4 = self.state_4.a
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
        self.a_5 = self.state_5.a

        self.state_5ss = H2(s = self.state_01.s,p = self.state_5.p)


        self.z_r = 1.5 * (self.b_5)
        self.n_r = np.round(np.pi/30*(110-self.alpha_4)*np.tan(self.alpha_4*np.pi/180)) #Glassman 1972
        self.q_5 = 2 * np.pi * self.r_5 / self.n_r
        self.o_5 = self.q_5 * self.C_m5 / self.W_5
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
        self.r_2 = 1.1 * self.r_3
        self.C_theta2 = self.C_theta1 * self.r_1/self.r_2
        C_m2_check = 0
        C_m2 = self.C_m1 * self.A_1 / self.A_2
        k=0
        while C_m2 != C_m2_check:
            if k !=0:
                C_m2 = C_m2_check
            # print("C_m2:",C_m2)
            self.C_2 = np.sqrt(C_m2**2 + self.C_theta2**2)
            self.h_2 = self.h_02 - 0.5 * self.C_2**2
            self.s_2 = self.state_1.s
            self.state_2 = H2(h=self.h_2,s=self.s_2)
            C_m2_check = self.CFE.mass_flow/(2 * np.pi * self.r_2 * self.state_2.rho * self.b_4)
            # print("C_m2 check:",C_m2_check)
            k+=1
        self.C_m2 = C_m2
        self.alpha_2 = np.tan(self.C_theta2/self.C_m2)*180/np.pi


        """Preliminary Calcs: Nozzle"""
        """Efficiency Correction Calulations"""
        self.h_loss = self.calc_losses()
        self.eta_ts_loss = self.deltah_0/(self.deltah_0 + self.h_loss)
        self.flow_coef = self.C_m5 / self.U_4
        self.load_coef = self.deltah_0 / (self.U_4**2)

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

    def print_turbine(self):
        turb_title = f'Turbine iteration: {self.i}\n\n'
        stn1_str_title = f'Station 1:\n\n'
        stn2_str_title = f'Station 2:\n\n'
        stn3_str_title = f'Station 3:\n\n'
        stn4_str_title = f'Station 4:\n\n'
        stn5_str_title = f'Station 5:\n\n'
        stn6_str_title = f'Station 6:\n\n'
        isen_vals_str = f' \
Total to static efficicency: {self.eta_ts}\n \
Stage loading: {self.load_coef}\n \
Flow coefficient: {self.flow_coef}\n \
Isentropic enthalpy drop: {self.h_0ss} [kJ kg^-1]\n \
Total-total enthalpy drop: {self.deltah_0} [kJ kg^-1]\n \
Total-to-static pressure ratio: {self.PR_ts}\n \
Specific speed: {self.N_s}\n \
Turbine velocity ratio: {self.v_s}\n \
Spouting Velocity: {self.C_0} [m s^-1]\n \
Turbine rotor inlet blade speed: {self.U_4} [m s^-1]\n\n'

        geom_str = f' \
Rotor inlet radius: {self.r_4*1000} [mm]\n \
Rotor inlet blade height: {self.b_4*1000} [mm]\n \
Rotor inlet relative flow angle: {self.beta_4} [degrees]\n \
Rotor outlet hub radius: {self.r_h5 * 1000} [mm]\n \
Rotor outlet shroud radius: {self.r_s5 * 1000} [mm]\n \
Rotor outlet blade height: {self.b_5 * 1000} [mm]\n \
Rotor axial length: {self.z_r * 1000} [mm]\n \
Rotor leading edge thickness: {self.t_lead * 1000} [mm]\n \
Rotor trailing edge thickness: {self.t_trail * 1000} [mm]\n \
Number of rotor blades: {self.n_r}\n \
Stator outlet radius: {self.r_3 * 1000} [mm]\n \
Stator inlet radius: {self.r_2 * 1e3} [mm]\n\n'
        C_theta4_opt = self.U_4 * (1-np.sqrt(np.cos(self.beta_4_rad))/(self.n_r**0.7))
        stn4_str = f' \
Rotor inlet stagnation pressure: {self.P_04/1000} [kPa]\n \
Rotor inlet radius: {self.r_4*1000} [mm]\n \
Absolute meridional velocity: {self.C_m4} [m s^-1]\n \
Absolute tangential velocity: {self.C_theta4} [m s^-1]\n \
Optimal absolute tangential velocity: {C_theta4_opt} [m s^-1]\n \
Absolute velocity: {self.C_4} [m s^-1]\n \
Rotor inlet absolute flow angle: {self.alpha_4} [degrees]\n \
Rotor inlet relative flow angle: {self.beta_4} [degrees]\n\n'
        print(turb_title)
        print(isen_vals_str + geom_str + stn4_str_title+stn4_str)

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
        

def find_turb(init_cfe,init_turb):
    new_turb = init_turb
    new_cfe = init_cfe
    i = new_turb.i
    while np.round(new_turb.eta_ts,4) != np.round(new_turb.eta_ts_loss,4):
        turb = new_turb
        cfe = new_cfe
        print("Turbine iteration:",turb.i)
        print("Turbine efficiency:",turb.eta_ts)
        print("Turbine specific speed:",turb.N_s)
        print("Turbine total-to-total enthalpy drop:",turb.deltah_0)
        print("Turbine total-to-static pressure ratio:",turb.PR_ts)
        print("Turbine enthalpy losses:",turb.h_loss)
        i += 1
        new_h_01 = (turb.eta_ts*turb.state_5ss.h - turb.h_05) / (turb.eta_ts - 1)
        print(new_h_01)
        new_state = H2(h=new_h_01,t=turb.T_01)
        new_P_01 =new_state.p
        new_PR = new_P_01 / turb.state_5.p

        new_dynamic_turb_inputs = {
            "PR_ts" : new_PR,
            "eta_ts" : turb.eta_ts_loss,
            "h_0ss" : 0,
        }

        new_cfe = CFE(static_cfe_inputs,new_dynamic_turb_inputs,i)
        new_delta_h0 = new_cfe.work_rate / new_cfe.mass_flow
        print(new_delta_h0)
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
        "PR_ts" : 1.002,
        "eta_ts" : 0.9,
        "h_0ss" : 0,
        "N_s" : 0
    }

    test_cfe = CFE(static_cfe_inputs,dynamic_turb_inputs,1)

    init_turb = turbine(test_cfe,test_cfe.static_turb_inputs,dynamic_turb_inputs,1)

    test_turb = find_turb(test_cfe,init_turb)
    test_turb.print_turbine()
    test_turb.velocity_triangles()
    # test_turb.print_states()
    # """"prelim calcs"""

    """Checks:"""
    print()
    state_5ss = H2(s = test_turb.state_01.s,p = test_turb.state_5.p)
    eff_check = (test_turb.h_01 - test_turb.h_05) / (test_turb.h_01 - state_5ss.h)
    print(test_turb.h_02 - test_turb.h_05)
    print(test_turb.deltah_0)
    print(eff_check)
    print()
    print("Mach No. Check:")
    M_inlet = test_turb.C_4/test_turb.state_4.a
    M_outlet = test_turb.C_5/test_turb.state_5.a
    print("Inlet Mach No.:",M_inlet)
    print("Outlet Mach No.",M_outlet)

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
