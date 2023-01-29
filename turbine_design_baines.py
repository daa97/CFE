import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.interpolate as inter
import csv
import time
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
        m_case = self.L * np.pi * (self.R_i**2 - .049**2) * rho_case

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

class station:
    def __init__(self,turb):
        pass
class whitfield_turbine:
    def __init__(self,CFE,static_turb_inputs,dynamic_turb_inputs,i):
        self.i = i
        self.state_01 = H2(t=static_turb_inputs["T_01"], p = self.CFE.P_in)
        self.T_01 = static_turb_inputs["T_01"]
        self.h_01 = self.state_01.h
        self.delta_W = self.CFE.work_rate / self.CFE.mass_flow
        self.S = self.delta_W / self.h_01 #Power ratio
        self.T_03 = (1 - self.S) * self.T_01

class turbine:
    def __init__(self,static_turb_inputs,dynamic_turb_inputs,i):
        self.CFE = CFE
        self.static_turb_inputs = static_turb_inputs
        self.i = i
        print(f'Turbine iteration {self.i}')
        #self.number_stations = static_turb_inputs["number of stations"]
        self.stations = self.make_stations()
        
        #self.eta_tt_guess = dynamic_turb_inputs["total to total efficiency guess"]
        self.eta_ts_guess = dynamic_turb_inputs["eta_ts"]
        self.PR_ts = dynamic_turb_inputs["PR_ts"]
        self.work_rate = static_turb_inputs["work_rate"]
        self.mass_flow = static_turb_inputs["mass_flow"]
        self.omega = static_turb_inputs["omega"]

        self.deltah_0 = self.work_rate/self.mass_flow
        self.T_01 = static_turb_inputs["T_01"]
        self.state_01 = H2(t=self.T_01, p = static_turb_inputs["P_01"])
        self.h_01 = self.state_01.h
        self.h_05 = self.state_01.h - self.deltah_0
        self.state_05 = H2(h = self.h_05, P=static_turb_inputs["P_01"]/self.PR_ts)
        self.epsilon_b = 0.4e-3
        self.Q = self.mass_flow / self.state_05.rho
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
            self.N_s = self.omega * np.sqrt(self.Q) / (self.h_0ss**(3/4))

        """Preliminary Calcs: Rotor"""

        """Isentropic Values"""
        self.v_s = self.calc_vel_ratio(dynamic_turb_inputs) # 0.737 * self.N_s**0.2
        self.C_0 = np.sqrt(2 * self.h_0ss)
        """Total and static states at station 4"""
        self.U_4 = self.C_0 * self.v_s # Blade spped [m/s]
        self.r_4 = self.U_4/self.omega # Rotor inlet radius [m]
        self.alpha_4 = 90 - 1* (10.8 + 14.2*self.N_s**2) # Rotor inlet angle [degrees]
        self.alpha_4_rad = self.alpha_4 / 180 * np.pi # Rotor inlet angle [radians]
        self.P_04 = self.state_01.p - self.state_01.rho * self.h_0ss*(1-self.eta_ts)/4 # Rotor inlet total pressure [Pa]
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
        self.b_4 = self.mass_flow / (2*np.pi*self.r_4*self.state_4.rho*self.C_m4)
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
        self.A_5 = self.mass_flow/(self.state_5.rho * self.C_m5)
        self.r_s5 = np.sqrt((self.A_5/np.pi)+self.r_h5**2) #Rotor outlet shroud radius
        self.r_5 = (self.r_s5+self.r_h5)/2 #Rotor outlet mean radius
        self.b_5 = self.r_s5-self.r_h5 #Rotor outlet blade height
        self.U_5 = self.omega * self.r_5
        self.W_m5 = self.C_m5
        self.W_theta5 = - self.U_5
        self.W_5 = np.sqrt(self.W_m5**2 + self.W_theta5**2)
        self.beta_5_rad = np.arctan(self.W_theta5/self.W_m5)
        self.beta_5 = self.beta_5_rad * 180 / np.pi
        self.P_05 = self.state_05.p
        self.state_5ss = H2(s = self.state_01.s,p = self.state_5.p)
        self.M_5 = self.C_5/ self.state_5.a
        self.M_5_rel = self.W_5 / self.state_5.a
        self.beta_s5_rad = np.arctan(self.r_5/self.r_s5 * np.tan(self.beta_5_rad))
        self.beta_s5 = self.beta_s5_rad * 180 / np.pi
        self.beta_h5_rad = np.arctan(self.r_5/self.r_h5 * np.tan(self.beta_5_rad))
        self.beta_h5 = self.beta_h5_rad * 180 / np.pi

        self.z_r = 1.5 * (self.b_5)
        self.n_r = np.round(np.pi/30*(110-self.alpha_4)*np.tan(self.alpha_4*np.pi/180)) #Glassman 1972
        self.q_5 = 2 * np.pi * self.r_5 / self.n_r
        self.o_5 = self.q_5 * self.C_m5 / self.W_5
        self.r_3 = self.r_4 + 2 * self.b_4 * np.cos(self.alpha_4*np.pi/180)
        self.r_2 = 1.1 * self.r_3
        self.t_lead = 0.04 * self.r_4
        self.t_trail = 0.02 * self.r_4
        """Preliminary Calcs: The Bowl"""
        self.r_1o = static_turb_inputs["R_o"]
        self.r_1i = self.r_1o - self.b_4
        self.r_cbo = self.r_1o - self.r_3
        self.r_cbi = self.r_1i - self.r_3
        # print(self.r_cbo)
        # print(self.r_cbi)

        """Efficiency Correction Calulations"""
        self.h_loss = self.calc_losses()
        self.eta_ts_loss = self.deltah_0/(self.deltah_0 + self.h_loss[0])
        self.flow_coef = self.C_m5 / self.U_4
        self.load_coef = self.deltah_0 / (self.U_4**2)
        self.D_s = 2 * self.r_4 * (self.h_0ss)**0.25 / self.Q**0.5
        rotor_eff = (self.state_04.h - self.state_05.h) / (self.state_04.h - self.state_5ss.h)
        
    def turbine_feasibility_checks(self):
        results = []
        VR1 = self.C_m5/self.U_4
        rad_ratio = self.r_h5/self.r_4
        rohlik_rad_ratio_lim = 1.29*self.N_s
        MVR = self.C_m5/self.C_m4
        R = (self.h_4 - self.h_5) / (self.h_01 - self.h_05)
        delta_zr_check = 1.5*self.b_4
        r_s5_check = 0.9 * self.r_4

        print(f'Exit meridional velocity ratio: {VR1}')
        if 0.2 <= VR1 and VR1 <= 0.4:
            print("Passed\n")
        else:
            print("Failed\n")

        print(f'Radius ratio: {rad_ratio}')
        if self.r_s5/self.r_4 <= 0.78:
            print("Baljie radius ratio test passed")
        else:
            print("Baljie radius ratio test failed")
        
        if rad_ratio <= 0.7:
            print("Rohlik radius ratio test passed")
        else:
            print("Rohlik radius ratio test failed")

        if self.r_s5 < r_s5_check:
            print("Aungier radius check passed\n")
        else:
            print("Aungier radius check failed\n")
        
        print(f"Meridional velocity ratio: {MVR}")
        if 1 <= MVR and MVR <= 1.5:
            print("Passed\n")
        else:
            print("Failed\n")
        
        print(f'Stage reaction: {R}')
        if 0.45 <= R and R<=0.65:
            print("Passed\n")
        else:
            print("Failed\n")
        
        print(f"Relative flow inlet angle: {self.beta_4}")
        if -20 >= self.beta_4 and self.beta_4 >= -40:
            print("Passed\n")
        else:
            print("Failed\n")
        
        print(f'Rotor axial length: {self.z_r}')
        if self.z_r >= delta_zr_check:
            print("Passed\n")
        else:
            print("Failed\n")

        blade_blockage_inlet = self.n_r * self.t_lead / (2 * np.pi * np.sin(np.pi/2 - self.beta_4_rad))

        print(f"Inlet blade metal blockage factor: {blade_blockage_inlet}")
        if blade_blockage_inlet < 0.5:
            print("Passed\n")
        else:
            print("Failed\n")        

        blade_blockage_outlet = self.n_r * self.t_trail / (2 * np.pi * np.sin(np.pi/2 - self.beta_5_rad))

        print(f"Outlet blade metal blockage factor: {blade_blockage_outlet}")
        if blade_blockage_outlet < 0.5:
            print("Passed\n")
        else:
            print("Failed\n")        

    def make_stations(self):#use this for making an h-s diagram later
        """Could I technically make all variables equal to -1 and then
        calculate them as their dependents become available?"""
        pass    

    def iterate_eta_ts(self):
        rho_05 = self.state_05.rho
        Q = self.mass_flow / rho_05
        eta_ts = 0
        """Eta Iteration"""
        while round(eta_ts,4) != round(self.eta_ts_guess,4):
            if eta_ts != 0:
                self.eta_ts_guess = eta_ts

            h_0ss = self.deltah_0 / self.eta_ts_guess
            
            N_s = self.omega * np.sqrt(Q) / (h_0ss**(3/4))

            eta_ts = 0.87 - 1.07 * (N_s - 0.55)**2 - 0.5 * (N_s - 0.55)**3
            # print("Calculated eta:", eta_ts)
            # print("Previous eta",self.eta_ts_guess)

        h_0ss = self.deltah_0 / eta_ts
        N_s = self.omega * np.sqrt(Q) / (h_0ss**(3/4))
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
        L_w = K_f * (rho_avg * self.U_4**3 * self.r_4**2) / (2 * self.mass_flow * self.W_5**2)
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
        # print(L_i/self.h_0ss,L_p/self.h_0ss,L_t/self.h_0ss,L_e/self.h_0ss,L_n/self.h_0ss,L_w/self.h_0ss)
        return [h_losses,L_e,L_i,L_p,L_w,L_t,L_n]

    def print_turbine(self,opts):
        turb_title = f'Turbine iteration: {self.i}\n\n'
        prelim_title = f'Preliminary Calculated Values:\n\n'
        dim_title = f'Non-dimensional Parameters:\n\n'
        geom_title = f'Turbine Geometry:\n\n'
        loss_title = f'Enthalpy Losses:\n\n'
        stn1_title = f'Station 1:\n\n'
        stn2_title = f'Station 2:\n\n'
        stn3_title = f'Station 3:\n\n'
        stn4_title = f'Station 4:\n\n'
        stn5_title = f'Station 5:\n\n'
        stn6_title = f'Station 6:\n\n'

        loss_str =f' \
Exit losses: {self.h_loss[1]/self.h_0ss}\n \
Incidence losses: {self.h_loss[2]/self.h_0ss}\n \
Passage losses: {self.h_loss[3]/self.h_0ss}\n \
Windage losses: {self.h_loss[4]/self.h_0ss}\n \
Trailing edge losses: {self.h_loss[5]/self.h_0ss}\n \
Nozzle losses: {self.h_loss[6]/self.h_0ss}\n \
Total losses: {self.h_loss[0]/self.h_0ss}\n \
    '
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
            "stations" : [stn4_title + stn4_str,stn5_title + stn5_str],
            "losses" : loss_title + loss_str
        }
        
        for key in opts:
            if key != "stations" and  (opts[key] == "y" or opts[key] == "Y"):
                print(strings[key])
            else:
                for i,stn in enumerate(opts[key]):
                    if stn =="Y" or stn == "y":
                        print(strings[key][i])

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

        """Calculate meridians"""
        m_h = np.zeros((N_p + N_p_L,1))
        hub = np.concatenate((m_h,hub_points),axis=1)
        m_s = np.zeros((N_p + N_p_L,1))
        shroud = np.concatenate((m_h,shroud_points),axis=1)

        for i in range(len(m_h)):
            if i==0:
                hub[i,0]=0
                shroud[i,0]=0
            else:
                hub[i,0] = np.sqrt((hub_points[i,0]-hub_points[i-1,0])**2 + (hub_points[i,1]-hub_points[i-1,1])**2) + hub[i-1,0]
                shroud[i,0] = np.sqrt((shroud_points[i,0]-shroud_points[i-1,0])**2 + (shroud_points[i,1]-shroud_points[i-1,1])**2) + shroud[i-1,0]    

        m_h = np.linspace(0,hub[-1,0],len(hub))
        m_s = np.linspace(0,shroud[-1,0],len(hub))
        hub_equal_coords = np.zeros((len(m_h),2))
        shroud_equal_coords = np.zeros((len(m_h),2))

        f_h = inter.interp1d(hub[:,0],hub[:,1:3],axis=0)
        f_s = inter.interp1d(shroud[:,0],shroud[:,1:3],axis=0)
        mid = np.zeros((len(m_h),3))


        for i,m in enumerate(m_h):
            if i==0:
                hub_equal_coords[i,0] = hub[0,1]
                hub_equal_coords[i,1] = hub[0,2]
                shroud_equal_coords[i,0] = shroud[0,1]
                shroud_equal_coords[i,1] = shroud[0,2]

                dx = (hub_equal_coords[i,0] - shroud_equal_coords[i,0])/2
                dy = (shroud_equal_coords[i,1] - hub_equal_coords[i,1])/2
                mid[i,1] = hub_equal_coords[i,0] - dx
                mid[i,2] = hub_equal_coords[i,1] + dy
                mid[i,0] = 0

            else:
                hub_equal_coords[i,0] = f_h(m_h[i])[0]
                hub_equal_coords[i,1] = f_h(m_h[i])[1]
                shroud_equal_coords[i,0] = f_s(m_s[i])[0]
                shroud_equal_coords[i,1] = f_s(m_s[i])[1]

                dx = (hub_equal_coords[i,0] - shroud_equal_coords[i,0])/2
                dy = (shroud_equal_coords[i,1] - hub_equal_coords[i,1])/2
                mid[i,1] = hub_equal_coords[i,0] - dx
                mid[i,2] = hub_equal_coords[i,1] + dy
                mid[i,0] = np.sqrt((mid[i,1]-mid[i-1,1])**2 + (mid[i,2]-mid[i-1,2])**2) + mid[i-1,0]

            plt.plot([hub_equal_coords[i,0],shroud_equal_coords[i,0]],[hub_equal_coords[i,1],shroud_equal_coords[i,1]])
        print(mid)
        plt.plot(mid[:,1],mid[:,2])
        plt.plot(shroud_points[:,0],shroud_points[:,1])
        plt.xlim([-.01, 0.06])
        plt.ylim([-.01, 0.06])
        plt.show()

        """Hub and shroud wrap angles"""
        beta_h5_aung = np.pi/2 - np.abs(self.beta_h5_rad)
        beta_s5_aung = np.pi/2 - np.abs(self.beta_s5_rad)
        beta_4_aung = np.pi/2
        cotbeta4 = 1/np.tan(beta_4_aung)
        cotbetas5 = 1/np.tan(beta_s5_aung)
        cotbetah5 = 1/np.tan(beta_h5_aung)
        print(beta_4_aung)
        print(beta_s5_aung)
        print(beta_h5_aung)
        print(cotbeta4)
        print(cotbetas5)
        print(cotbetah5)

        m_4 = mid[-1,0]

        theta_h = np.zeros((len(hub_points),))
        theta_s = np.zeros((len(hub_points),))
        norm_m = np.zeros((len(hub_points),))
        norm_m_h = np.zeros((len(hub_points),))
        m = mid[:,0]
        # cot_beta4 = 1/np.tan()
        theta_4 = (m_4 / 2) * (cotbeta4 / self.r_4 + cotbetas5 / self.r_s5)
        print(theta_4*180/np.pi)

        A = cotbetas5 / self.r_s5

        B = m_4**(-2) * (cotbeta4/self.r_4 - A)

        C = - B / (2 * m_4)

        D = cotbetah5 / self.r_h5

        E = (3 * theta_4/ m_4**2) - 1/m_4 * (2 * D + cotbeta4 / self.r_4)

        F = m_4**(-2) * (D + cotbeta4 / self.r_4) - 2 * theta_4 / m_4**3
        print(A,B,C,D,E,F)
        for i,theta in enumerate(theta_h):

            theta_s[i] = (A * m[i] + B * m[i]**3 + C * m[i]**4) 
            # print("shroud",m_s[i])
            # print("hub",m_h[i])
            theta_h[i] = (D * m[i] + E * m[i]**2 + F * m[i]**3) 
            norm_m[i] = mid[i,0]/m_4

        beta_h = np.zeros(len(theta_h),)
        print(beta_h)
        beta_s = np.zeros(len(theta_h),)
        for i,beta in enumerate(beta_h):
            if i==0:
                beta_h[i] = beta_h5_aung * 180/np.pi
                beta_s[i] = beta_s5_aung * 180/np.pi
            if i==1:
                pass
            else: 
                beta_h[i-1] = hub[i-1,2] * (theta_h[i] - theta_h[i-2]) / (mid[i,0]-mid[i-2,0]) * 180/np.pi
                beta_s[i-1] = shroud[i-1,2] * (theta_s[i] - theta_s[i-2]) / (mid[i,0]-mid[i-2,0]) * 180/np.pi
        beta_h[-1] = beta_4_aung * 180/np.pi
        beta_s[-1] = beta_4_aung * 180/np.pi
        # print("hub meridion")
        # print(m_h)
        # print("shroud meridion")
        # print(m_s)
        # print("Hub wrap angle")
        # print(theta_h)
        # print("Shroud wrap angle")
        # print(theta_s)
        plt.plot(norm_m,theta_h,marker=".")
        plt.plot(norm_m,theta_s,marker=".")
        plt.show()
        plt.plot(norm_m,beta_h)
        plt.plot(norm_m,beta_s)
        plt.show()

    def calc_vel_ratio(self,dynamic_turb_inputs):
        if dynamic_turb_inputs["v_s"] == "default":
            return 0.737 * self.N_s**0.2
        else:
            return dynamic_turb_inputs["v_s"]

class nozzle:
    def __init__(self,nozzle_inputs,turb):
        self.N_n = 13 # Number of nozzle blades guess
        self.N_p_n = 20
        self.turb = turb
        self.sc = nozzle_inputs["sc"]
        self.theta = nozzle_inputs["camber angle"]
        self.ac = nozzle_inputs["ac"]
        self.t_2c = nozzle_inputs["t_2c"]
        self.t_3c = nozzle_inputs["t_3c"]
        self.t_maxc = nozzle_inputs["t_maxc"]
        self.dc = nozzle_inputs["dc"]
        self.radius_ratio = nozzle_inputs["radius ratio"]
        self.gam_guess = nozzle_inputs["setting angle"]
        # self.chi_2 = np.arctan((4 * self.bc)/(4*self.ac -1))
        # self.chi_3 = np.arctan()
        self.r_3 = turb.r_3
        self.r_2 = self.r_3 * self.radius_ratio
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
        self.alpha_3_rad =np.pi/2 - np.abs(np.arctan(self.C_theta3/self.C_m3))
        self.alpha_3 = self.alpha_3_rad * 180 / np.pi
        print(self.alpha_3)
        self.d_3 = 2 * np.pi * self.r_3 / self.N_n # Pitch is usually s but since s is entropy I am using d
        print(self.d_3)
        self.o_3 = self.d_3 * np.sin(self.alpha_3_rad)
        print(self.o_3)
        self.c = 1/self.sc * self.d_3

        naca_output = self.calc_naca_profile()
        self.xcs = naca_output[0]
        self.ycs = naca_output[1]
        self.naca_cam = np.array([self.xcs,self.ycs],dtype=object).T
        self.naca_chi = naca_output[2]
        self.naca_suc_surf = naca_output[3]
        self.naca_pres_surf = naca_output[4]

        self.gamma_3 = self.find_setting_angle()
        self.gamma_2 = np.arccos(self.r_3 * np.cos(self.gamma_3) / self.r_2)
        print("Gamma 3:",self.gamma_3*180/np.pi)
        self.beta_3b = self.gamma_3 - self.naca_chi[-1]
        self.beta_2b = self.gamma_2 - self.naca_chi[0]
        print(self.beta_2b*180/np.pi)

        self.alpha_2 = self.find_alpha_2()
        
    def calc_naca_profile(self):
        N_p_n = self.N_p_n
        check = np.zeros((N_p_n,))
        tcs = np.zeros((N_p_n,))
        suc = np.zeros((N_p_n,2))
        pres = np.zeros((N_p_n,2))

        xcs = np.linspace(0,self.c,N_p_n)
        ycs = np.zeros((N_p_n,))
        print(np.tan(self.theta))
        print((self.ac - self.ac**2 - 3/16))
        print(1 + (4 * np.tan(self.theta))**2 * (self.ac - self.ac**2 - 3/16))
        print((4 * np.tan(self.theta)))
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
        # print(chi)
        return [xcs,ycs,chi,suc,pres]

    def create_cascade(self,gam_guess):
        rot_cam = np.zeros((self.N_p_n,2))
        rot_suc = np.zeros((self.N_p_n,2))
        rot_pres = np.zeros((self.N_p_n,2))
        blade_2_cam = np.zeros((self.N_p_n,2))
        blade_2_suc = np.zeros((self.N_p_n,2))
        blade_2_pres = np.zeros((self.N_p_n,2))
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
            if np.abs(blade_2_cam[i,0]) < np.abs(blade_2_pres[i,0]):
                blade_2_pres[i,0] = - blade_2_pres[i,0]
            blade_2_pres[i,1] = blade_2_pressure[1]

            blade_2_suc[i,0] = blade_2_suction[0]
            if np.abs(blade_2_cam[i,0]) > np.abs(blade_2_suc[i,0]):
                blade_2_suc[i,0] = - blade_2_suc[i,0]
            blade_2_suc[i,1] = blade_2_suction[1]


        """Find the throat of the stator cascade"""
        min_point = rot_suc[-1]
        d = np.zeros((self.N_p_n,))
        for i,point in enumerate(blade_2_pres):
            d[i] = np.sqrt((min_point[0]-point[0])**2 + (min_point[1]-point[1])**2)
            if i > 0:
                if d[i] < d[i-1]:
                    min_d = d[i]
                    min_d_point = point
        throatx = [min_point[0],min_d_point[0]]
        throaty = [min_point[1],min_d_point[1]]

        thetas = np.linspace(0, 2 * np.pi, num = 100)
        circ = np.zeros((100,2))
        for i,theta in enumerate(thetas):
            circ[i,0] = self.r_3 * (np.cos(theta))
            circ[i,1] = self.r_3 * (np.sin(theta))
        plt.plot(throatx,throaty,color="r")
        plt.plot(circ[:,0],circ[:,1])
        plt.plot(blade_2_suc[:,0],blade_2_suc[:,1])
        plt.plot(blade_2_pres[:,0],blade_2_pres[:,1])
        plt.plot(blade_2_cam[:,0],blade_2_cam[:,1])
        plt.plot(rot_cam[:,0],rot_cam[:,1])
        plt.plot(rot_suc[:,0],rot_suc[:,1])
        plt.plot(rot_pres[:,0],rot_pres[:,1])
        plt.show()

        return [rot_suc,rot_pres,rot_cam,min_d,min_point,min_d_point]

    def find_setting_angle(self):
        throat_width = self.o_3
        print("Throat width:",throat_width)
        gam = self.gam_guess
        print("gam:",gam)
        rotated_blade = self.create_cascade(gam)
        throat_width_guess = rotated_blade[3]
        print("o:",throat_width_guess)
        while np.round(throat_width,4) != np.round(throat_width_guess,4):
            gam = np.arcsin(np.sin(gam)*(throat_width/throat_width_guess))
            print("gamma:",gam)
            rotated_blade = self.create_cascade(gam)
            throat_width_guess = rotated_blade[3]
            print("o:",throat_width_guess)
        return gam

    def find_alpha_2(self):
        r = np.linspace(self.r_2,self.r_3,self.N_p_n)
        gammas = np.zeros((self.N_p_n,))
        betas = np.zeros((self.N_p_n,))
        gammas[0] = self.gamma_2
        betas[0] = self.beta_2b
        L = np.zeros((self.N_p_n))
        for i,gamma in enumerate(gammas):
            if i > 0:
                gammas[i] = np.arccos(self.r_3 * np.cos(self.gamma_3) / r[i])
                betas[i] = gammas[i] - self.naca_chi[i]
        print("Gammas")
        print(gammas)
        print("Betas")
        print(betas)
        for i in range(len(L)):
            if i > 0:
                L[i] = L[i-1] + (r[i]-r[i-1])/(np.sin((betas[i])))
        print(L)
        i_star_1 = (3.6 * np.sqrt(10*self.t_2c*self.c/L[-1]*-1) + np.abs((self.beta_3b-self.beta_2b)/3.4))
        i_star_2 = np.sqrt(-1*L[-1]/(self.sc*self.c))
        i_star_3 = np.abs(self.beta_3b-self.beta_2b)
        i_star = i_star_1*i_star_2-i_star_3
        print(i_star*180/np.pi)

        alpha_2 = self.beta_2b - i_star * np.sign(self.beta_3b-self.beta_2b)
        print(alpha_2*180/np.pi)

        return alpha_2

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
            "h_0ss" : new_delta_h0ss,
            "v_s" : turb.v_s
            }
        new_static_turb_inputs = turb.static_turb_inputs
        new_turb = turbine(new_static_turb_inputs,new_dynamic_turb_inputs,i)
    new_turb.turbine_feasibility_checks()
    return new_turb

def find_stator(nozzle_inputs,turb):
    # gams = np.zeros(self.N_pb)
    pass
    
if __name__ == "__main__":
    pass

