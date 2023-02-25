import numpy as np

class CFE:
    def __init__(self,**kwargs):
        self.static_cfe_inputs = kwargs
        """ Entry channel fluid properties """
        self.T_in = kwargs["temp"] #Temperature [K]
        self.P_in = kwargs["press"]
        self.fluid = kwargs["fluid"]
        self.cfe_state = self.fluid(p=self.P_in,t = self.T_in)
        self.gamma = self.cfe_state.gamma
        T = self.T_in
        self.mu = self.cfe_state.mu
        self.rho = self.cfe_state.rho #Density, [kg/m^3]
        self.nu = self.mu/self.rho #Kinematic viscosity [m^2 s^-1]
        
        """ Entry channel geometric properties """
        self.R_i = kwargs["inner_radius"] #CFE channel inner radius [m]
        self.R_o = kwargs["outer_radius"] #CFE channel outer radius [m]
        self.uranium_mass = kwargs["uranium_mass"]
        self.annulus_area = np.pi * (self.R_o**2 - self.R_i**2) #CFE channel annulus area [m^2]
        self.h = self.R_o - self.R_i #CFE channel thickness [m]
        self.eta = self.R_i/self.R_o #CFE channel radius ratio
        self.L = kwargs["length"] #CFE Channel length [m]
        
        """ CFE Turbine requirements and inputs"""
        self.mass_flow = kwargs["mass_flow"] #Mass Flow Rate [kg/s]
        self.omega = kwargs["rpm"]*np.pi/30 #Angular Velocity [s^-1]
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
            "omega" : self.omega,
            "fluid" : self.fluid
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
        A = 3.33602; B = -35.00831; C = 107.19577; # experimental bearing friction coefficients
        I=0.00039205; d=.06      # inertia and diameter of tested bearing
        base_load = 2*I/(fric_coeff * d/2) *np.sqrt(B**2 - 3*A*(C-self.omega))
        print("BASE:", self.omega, base_load)
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
        print(f"Re w, a: {Re_w:,.1f}, {Re_a:,.1f}")
        lR = np.log10(Re_w)
        lewis_eta = 0.15999/0.22085
        lewis_G_factor = 4*np.pi*lewis_eta/((1-lewis_eta)*(1-lewis_eta**2))
        C0 = np.log10(lewis_G_factor)
        if 250<=Re_w<=13e3:
            exp = .2005*lR**3 - 1.970*lR**2 + (7.775-1)*lR - 5.516 - C0
        elif 13e3<Re_w<=1e6:
            exp = -.006360*lR**3 + .1349*lR**2 + (.8850-1)*lR + 1.610 - C0
        else:
            raise ValueError("Rotational Reynolds number out of range!")
        lewis_Nu = 10**exp
        Nu = 1.1* lewis_Nu      # adjustment for axial flow
        #print(f"Nu_w: {Nu:.4f}")
        Mlam = 4*np.pi*self.cfe_state.mu*self.L*self.omega / (self.R_i**(-2) - self.R_o**(-2))
        #print(f"Wlam: {Mlam*self.omega:.4f}")
        self.M_visc = Nu*Mlam
        #print(f"W_visc: {self.M_visc*self.omega}")
        return self.M_visc

    def calc_inlet_inertial_moment(self):
        R_mean = (self.R_i+self.R_o) / 2
        Isp_fluid = R_mean**2           # specific moment of inertia of the fluid
        omega_ttcf = self.omega/2       # approx. fluid rotation rate at stator inlet
        pipe_rad = 0.493/2 * 25.4/1e3   # inlet manifold pipe radius
        pipe_area = np.pi * pipe_rad**2 # 
        vdot = self.mass_flow/self.rho
        #print("Vdot:", vdot)
        U_inlet = vdot / pipe_area
        #print("U inlet:", U_inlet)
        omega_inlet = U_inlet / R_mean
        #print("omega inlet:", omega_inlet)
        dL_sp = Isp_fluid * (omega_ttcf - omega_inlet) # change in angular momentum per unit mass of fluid
        #print("omega ttcf", omega_ttcf)
        #print("speed of sound:", self.cfe_state.A)
        #print("Mach No:", U_inlet/self.cfe_state.A)
        self.M_inlet = dL_sp * self.mass_flow
        if self.M_inlet < 0:
            self.M_inlet /= 2
        return self.M_inlet

    def calc_work_rate(self):
        """Calculates the required work rate of the turbine based on viscous
        losses from shear at the CFE surface. This calculation is an 
        amalgamation of various papers of torque from taylor-couette flow.
        The correlation for non-dimensional torque comes from Lewis 1999.
        """
        M_bearings = self.calc_bearing_moment()
        M_visc = self.calc_visc_moment()       
        M_inlet = self.calc_inlet_inertial_moment()
        M = M_bearings + M_visc + M_inlet
        print(f'''Bearing Moment: {M_bearings}
                Viscous Moment: {M_visc}
                Inlet Moment: {M_inlet}''')
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


