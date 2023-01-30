import numpy as np
np.set_printoptions(precision=1)
import matplotlib.pyplot as plt
import scipy as sp
import sys
from fluids import H2
import multiprocessing as multi


class Face:
    def __init__(self,face_radius, CFE_len):
        self.radius     = face_radius #m
        self.area       = 2 * np.pi * face_radius * CFE_len #m^2

class Cell:
    def __init__(self,cell_len,cell_temp,cell_radius,CFE_len,BCs,i,num_cells,press,m_b,mass_flow,omega): 
        self.press  = press
        self.cell_temp   = float(cell_temp) #FIXME - assign temp from temp profile
        self.H2     = H2(P=self.press, T=self.cell_temp, linear=True)
        self.rho_H  = self.H2.rho
        self.Cp     = self.H2.Cp
        self.h      = self.H2.h

        self.i           = i
        self.cell_len    = cell_len 
        self.BC          = BCs 
        self.cell_radius = cell_radius #Taken from CFE create_mesh method


        self.omega  = omega #rpm
        self.m_dot  = mass_flow #kg/s

        self.k_U    = 13.7 #W/m-K 
        self.rho_U  = self.calc_rho_U() #kg/m^3
        self.Cp_U   = 179.78 + 0.0131 * self.cell_temp #J/kg/K - R**2 - 0.9998


        self.V_b_o  = 1e-9 #m^3
        self.m_b    = m_b #kg
        self.m_b_o  = self.calc_init_bubble_mass() #kg
        self.V_B    = self.m_b_o * self.H2.volume #m^3
        self.u_B    = self.calc_avg_velocity() #m/s
        

        self.num_cells  = num_cells
        self.CFE_len    = CFE_len
        self.g          = self.calc_g()
        
        self.face_e = Face(self.cell_radius-self.cell_len/2,self.CFE_len)
        self.face_w = Face(self.cell_radius+self.cell_len/2,self.CFE_len)

        self.cell_vol   = np.pi * ((self.face_w.radius)**2 - (self.face_e.radius)**2) * self.CFE_len
        self.gas_vol    = self.m_dot * self.H2.volume * self.cell_len / self.u_B
        self.void       = self.calc_void_fraction()
        self.vol_void   = self.gas_vol / self.cell_vol

        self.num_b      = self.calc_num_b()
        self.drag       = self.calc_drag()

    def calc_avg_velocity(self):
        a_c = self.cell_radius * (self.omega ** 2)
        a = ((3 * self.V_B )/ (4 * np.pi)) ** (1/3)
        r_b = 2.179 * a
        u_B = 2/3 * np.sqrt(r_b * a_c)
        return u_B

    
    def calc_num_b(self):
        num_b = (self.m_dot/self.m_b_o) / self.u_B * self.cell_len
        return num_b

    def calc_g(self):
        g = self.omega ** 2 * self.cell_radius 
        return g

    def calc_void_fraction(self):
        N_B_dot = self.m_dot / self.cell_radius / self.u_B / self.rho_H
        void = N_B_dot / 2 / np.pi / self.CFE_len
        return void

    def calc_init_bubble_mass(self):
        if self.i ==0:
            m_b_o = self.rho_H * self.V_b_o
            return m_b_o
        else:
            return self.m_b

    def find_phase_U(self):
        if self.temp > 4300:
            phase = 'v'
        else:
            phase = 'l'
        return phase

    def calc_rho_U(self):
        rho = (19022 - 1.4486 * self.cell_temp)
        return rho

    def calc_drag(self):
        a = ((self.V_B * 3)/(4*np.pi))**(1/3)
        SA = np.pi * (2.179 * a * np.sin(52/180 * np.pi))**2
        drag_per = 4/3 * self.rho_U * self.g * a * SA/self.V_B * self.void
        drag = drag_per * self.cell_vol
        return drag


class CFE:
    def __init__(self,IR,annulus_t,length, power,num_cells,temp_profile,mass_flow,rpm,P_0,BC,q_tp,iteration=0,OG_power=0,press_profile=None):
        self.it = iteration
        self.P_0 = P_0
        self.annulus_t = annulus_t
        self.OR_U = IR #Uranium annulus outer radius in cm
        self.IR_U = (IR-annulus_t) #Uranium annulus inner radius
        self.OR_H = self.IR_U #Hydrogen outer radius
        self.L = length #Length of the CFE
        self.MW = power #Total CFE power in megawatts
        self.OG_MW = OG_power
        self.cell_len = annulus_t/num_cells
        self.mass_flow = mass_flow
        self.num_cells = num_cells
        self.cell_center_radii = self.create_cell_center_radii()
        self.temp_profile_input = temp_profile
        self.temp_profile = self.create_temp_profile()
        self.BC = BC
        self.rpm = rpm
        if self.it==0:
            self.press_profile = [self.P_0] * self.num_cells
        else:
            self.press_profile = press_profile
        self.mesh = self.create_mesh()
        self.q_profile_input = q_tp
        self.q_tp = self.create_heat_profile()
        self.sys = self.create_sys()
        self.matrix = self.sys[0]
        self.b = self.sys[1]
        self.q_tot = self.sys[2]
        self.X = self.sys[3]



    def create_cell_center_radii(self):
        cell_center_radii = []
        for i in range(self.num_cells):
            cell_center_radii.append(self.OR_U-self.cell_len/2-self.cell_len*i)
        return cell_center_radii

    def create_temp_profile(self):
        if type(self.temp_profile_input) == str:
            temp_profile = []
            for i in range(self.num_cells):
                r = self.cell_center_radii[i]
                temp_profile.append(eval(self.temp_profile_input))

        elif type(self.temp_profile_input) is list:
            temp_profile = self.temp_profile_input

        elif type(self.temp_profile_input) is np.ndarray:
            temp_profile = []
            for i in range(self.num_cells):
                temp_profile.append(float(self.temp_profile_input[i]))

        return temp_profile
    
        

    def create_heat_profile(self):
        if type(self.q_profile_input) == str:
            q_profile = []
            
            for i in range(self.num_cells):
                r = self.cell_center_radii[i]
                q_profile.append(eval(self.q_profile_input))
            
            norm = sum(q_profile)
            q_tp = [self.MW*(10**6)*q/norm for q in q_profile]

        elif type(self.q_profile_input) is list:
            q_tp = [self.MW*(10**6)*percent for percent in self.q_profile_input]

        elif type(self.q_profile_input is np.ndarray):
            q_t = np.polynomial.polynomial.polyval(np.multiply(100,self.cell_center_radii),np.multiply(self.MW/10,self.q_profile_input))
            cell_vol = [(100**3)*cell.cell_vol for cell in self.mesh]
            q_tp = np.multiply(q_t,cell_vol)
        return q_tp

    def find_press_profile(self):
        preint = (self.rpm * 2*np.pi/60)**2
        dPs = []
        Ps = []
        rho_total = lambda cell: cell.rho_U * (1-cell.void) + cell.rho_H * cell.void
        for c in list(reversed(self.mesh)):
            dP = c.cell_radius * rho_total(c) * c.cell_len * preint
            dPs.append(dP)
            Ps.append(np.sum(dPs)+self.P_0)
        
        return list(reversed(Ps))

    def create_mesh(self):
        mesh = []
        m_b = -1
        for i in range(self.num_cells):
            if i == 0:
                BCs = self.BC[0]
            elif i == self.num_cells-1:
                BCs = self.BC[1]
            else:
                BCs = [0,"none"]
            mesh.append(Cell(cell_len = self.cell_len,cell_temp = self.temp_profile[i],cell_radius=self.cell_center_radii[i],CFE_len = self.L,
                        BCs = BCs,i=i,num_cells=self.num_cells,press=self.press_profile[i],m_b = m_b,mass_flow = self.mass_flow,omega = self.rpm*2 *np.pi/60))
            if i == 0:
                cell_0 = Cell(cell_len = self.cell_len,cell_temp = self.temp_profile[i],cell_radius=self.cell_center_radii[i],CFE_len = self.L,BCs = BCs,i=i,num_cells=self.num_cells,press=self.press_profile[i],m_b = m_b,mass_flow = self.mass_flow,omega = self.rpm*2 *np.pi/60)
                m_b = cell_0.m_b_o
        return mesh

    def create_sys(self): #Since Pe < 2 uses central differencing for discretization of convection term and a pseudo-transient for stability
        dt = 0.01
        q_tot = 0
        sys = np.zeros((self.num_cells,self.num_cells),float)
        b = np.zeros((self.num_cells,1))
        X=[]
        for i,cell in enumerate(self.mesh):
            X.append(cell.void)

            a_c =(1-cell.void) * cell.k_U * (cell.face_w.area + cell.face_e.area) / cell.cell_len +(1-cell.void) * cell.rho_U * cell.Cp_U * cell.cell_vol / dt
            try:
                a_w = - (1-cell.void) * cell.k_U * cell.face_w.area / cell.cell_len
            except IndexError:
                a_e = - (1-cell.void) * cell.k_U * cell.face_e.area / cell.cell_len 

            try:
                a_e =  - (1-cell.void) * cell.k_U * cell.face_e.area / cell.cell_len 
            except IndexError:
                a_w = - (1-cell.void) * cell.k_U * cell.face_w.area / cell.cell_len

            SF_v = (1 - cell.void)/.5
            q_tot += SF_v * self.q_tp[i]

            if i == 0 or i == (self.num_cells - 1):
                if cell.BC[1].lower() == 'dirichlet':
                    if i == 0:
                        T_ghost =  2 * cell.BC[0] - cell.cell_temp #FIXME - Ghost Dirichlet BC is causing converging oscillations
                        ghost = Cell(self.cell_len,T_ghost,cell.cell_radius+cell.cell_len,self.L,BCs=[],i='',num_cells=self.num_cells,press=self.P_0,m_b = cell.m_b_o,mass_flow = self.mass_flow,omega = self.rpm*2 *np.pi/60)
                        a_w_d =  - (1-cell.void) * cell.k_U * cell.face_w.area / cell.cell_len
                        sys[i][i] = a_c - a_w_d
                        sys[i][i+1] = a_e
                        b[i] = self.q_tp[i] * SF_v - 2 * a_w_d * cell.BC[0] - self.mass_flow/2 * (self.mesh[i+1].h - ghost.h) + (1-cell.void) * cell.rho_U * cell.Cp_U * cell.cell_vol * cell.cell_temp / dt
                    else:
                        pass
                elif cell.BC[1].lower() == 'neumann':
                    if i == 0:
                        pass
                    else:
                        a_e_n = - (1-cell.void) * cell.k_U * cell.face_e.area / cell.cell_len
                        sys[i][i] = a_e_n + a_c
                        sys[i][i-1] = a_w
                        b[i] = self.q_tp[i] * SF_v - a_e_n * cell.BC[0] * cell.cell_len - self.mass_flow/2 * (cell.h - self.mesh[i-1].h) + (1-cell.void) * cell.rho_U * cell.Cp_U * cell.cell_vol * cell.cell_temp / dt
            else:
                sys[i][i] = a_c
                sys[i][i-1] = a_w
                sys[i][i+1] = a_e
                b[i] = self.q_tp[i] * SF_v - self.mass_flow/2 * (self.mesh[i+1].h - self.mesh[i-1].h) + (1-cell.void) * cell.rho_U * cell.Cp_U * cell.cell_vol * cell.cell_temp / dt

        if np.abs(self.OG_MW*(10**6) - q_tot) >1e-4:
            mod = (self.OG_MW*(10**6))/q_tot
            self.MW = self.MW * mod

        return [sys,b,q_tot,X]

def status(msg):
    sys.stdout.write("\r"+msg+" "*20)
    sys.stdout.flush()

def solve(init_CFE):
    res = 100
    maxres = res
    tol = 1e-6
    c = 0
    old_CFE = init_CFE
    T_old = np.array(old_CFE.temp_profile).reshape(len(old_CFE.temp_profile),1)
    while  res > tol:
        if maxres == 100:           # get initial residual value
            maxres = res
        #status(f"Err: {res/tol:.2f}; Completion: {(1-np.log(res/tol)/np.log(maxres/tol))*100:.3f}%")    # print progress 

        next_CFE = CFE(old_CFE.OR_U, old_CFE.annulus_t, old_CFE.L, old_CFE.MW, old_CFE.num_cells, T_old,old_CFE.mass_flow,
                        old_CFE.rpm, old_CFE.P_0,old_CFE.BC,old_CFE.q_profile_input,iteration = old_CFE.it + 1, OG_power = old_CFE.OG_MW, press_profile=old_CFE.find_press_profile())
        T_new = sp.linalg.solve(next_CFE.matrix,next_CFE.b)
        res = np.amax(np.divide(np.absolute(T_new-T_old),T_old))    

        T_old = T_new[:]
        old_CFE = next_CFE
        c += 1
        if c > 10000:
            res = 0
            
    next_CFE = CFE(old_CFE.OR_U,old_CFE.annulus_t, old_CFE.L, old_CFE.MW, old_CFE.num_cells, T_old,old_CFE.mass_flow, old_CFE.rpm, old_CFE.P_0, old_CFE.BC, old_CFE.q_profile_input,iteration = old_CFE.it + 1,OG_power=old_CFE.OG_MW, press_profile=old_CFE.find_press_profile())
    xi = np.array([next_CFE.BC[0][0], next_CFE.P_0])
    wall_h = H2(T=xi[0],P=xi[1], linear=True).h
    e_to_H2 = next_CFE.mass_flow * (next_CFE.mesh[next_CFE.num_cells-1].h - wall_h)

    T_ghost = 2 * next_CFE.mesh[0].BC[0] - next_CFE.mesh[0].cell_temp #FIXME - Ghost Dirichlet BC is causing converging oscillations
    ghost = Cell(next_CFE.cell_len,T_ghost,next_CFE.mesh[0].cell_radius+next_CFE.mesh[0].cell_len,next_CFE.L,BCs=[],i='',num_cells=next_CFE.num_cells,press=next_CFE.P_0,m_b = -1,mass_flow = next_CFE.mass_flow,omega = next_CFE.rpm*2 *np.pi/60)    
    mom=[]
    fld_mom = []
    fld_arch=[]
    fld_drag=[]
    fld_cont=[]
    fld_u = []
    rho = []
    rho_H = []
    for i,cell in enumerate(next_CFE.mesh):
        fld_u.append(cell.u_B)
        rho.append(cell.rho_U)
        rho_H.append(cell.rho_H)
        if i == 0:
            cell_mom = (next_CFE.mass_flow / 2  * (next_CFE.mesh[i+1].void * next_CFE.mesh[i+1].u_B - ghost.void * ghost.u_B))
            fld_mom.append(cell_mom)
            cell_arch = (cell.void * cell.cell_vol * cell.g * (cell.rho_H-cell.rho_U)) 
            fld_arch.append(cell_arch)
            fld_drag.append(cell.drag)
            fld_cont.append(1 / 2  * (next_CFE.mesh[i+1].void * next_CFE.mesh[i+1].rho_H * next_CFE.mesh[i+1].u_B * cell.face_e.area   - ghost.void * ghost.rho_H * ghost.u_B * cell.face_w.area))
            add = (next_CFE.mass_flow / 2  * (next_CFE.mesh[i+1].void * next_CFE.mesh[i+1].u_B - ghost.void * ghost.u_B)) + (cell.void * cell.cell_vol * cell.g * (cell.rho_H-cell.rho_U)) + cell.drag
            mom.append(add)

        elif i == 1:
            cell_mom = next_CFE.mass_flow * ((3/2 * cell.void * cell.u_B) - (1/2 * next_CFE.mesh[i-1].void * next_CFE.mesh[i-1].u_B) + ((-3/2 * next_CFE.mesh[i-1].void) + (1/2 * ghost.void * ghost.u_B)))
            fld_mom.append(cell_mom)
            cell_arch = (cell.void * cell.cell_vol * cell.g * (cell.rho_H-cell.rho_U)) 
            fld_arch.append(cell_arch)
            fld_drag.append(cell.drag)
            fld_cont.append(1/2 * (next_CFE.mesh[i+1].void * next_CFE.mesh[i+1].rho_H * next_CFE.mesh[i+1].u_B * cell.face_e.area - next_CFE.mesh[i-1].void * next_CFE.mesh[i-1].rho_H * next_CFE.mesh[i-1].u_B*cell.face_w.area))
            mom.append(next_CFE.mass_flow * ((3/2 * cell.void * cell.u_B) - (1/2 * next_CFE.mesh[i-1].void * next_CFE.mesh[i-1].u_B) + ((-3/2 * next_CFE.mesh[i-1].void) + (1/2 * ghost.void * ghost.u_B))) + (cell.void * cell.cell_vol * cell.g * (cell.rho_H-cell.rho_U))+ cell.drag)
        
        elif i == next_CFE.num_cells-1:
            cell_mom = next_CFE.mass_flow * ((3/2 * cell.void * cell.u_B) - (1/2 * next_CFE.mesh[i-1].void * next_CFE.mesh[i-1].u_B) + ((-3/2 * next_CFE.mesh[i-1].void)+(1/2 * next_CFE.mesh[i-2].void * next_CFE.mesh[i-2].u_B))) 
            fld_mom.append(cell_mom)
            cell_arch = (cell.void * cell.cell_vol * cell.g * (cell.rho_H-cell.rho_U)) 
            fld_arch.append(cell_arch)
            fld_drag.append(cell.drag)
            fld_cont.append(1/2 * (next_CFE.mesh[i].void * next_CFE.mesh[i].rho_H * next_CFE.mesh[i].u_B * cell.face_e.area - next_CFE.mesh[i-1].void * next_CFE.mesh[i-1].rho_H * next_CFE.mesh[i-1].u_B*cell.face_w.area))
            mom.append(next_CFE.mass_flow * ((3/2 * cell.void * cell.u_B) - (1/2 * next_CFE.mesh[i-1].void * next_CFE.mesh[i-1].u_B) + ((-3/2 * next_CFE.mesh[i-1].void)+(1/2 * next_CFE.mesh[i-2].void * next_CFE.mesh[i-2].u_B))) + (cell.void * cell.cell_vol * cell.g * (cell.rho_H-cell.rho_U)) + cell.drag)
   
        else:
            cell_mom = next_CFE.mass_flow * ((3/2 * cell.void * cell.u_B) - (1/2 * next_CFE.mesh[i-1].void * next_CFE.mesh[i-1].u_B) + ((-3/2 * next_CFE.mesh[i-1].void)+(1/2 * next_CFE.mesh[i-2].void * next_CFE.mesh[i-2].u_B))) 
            fld_mom.append(cell_mom)
            cell_arch = (cell.void * cell.cell_vol * cell.g * (cell.rho_H-cell.rho_U)) 
            fld_arch.append(cell_arch)
            fld_drag.append(cell.drag)
            fld_cont.append(1/2 * (next_CFE.mesh[i+1].void * next_CFE.mesh[i+1].rho_H * next_CFE.mesh[i+1].u_B * cell.face_e.area - next_CFE.mesh[i-1].void * next_CFE.mesh[i-1].rho_H * next_CFE.mesh[i-1].u_B*cell.face_w.area))

            mom.append(next_CFE.mass_flow * ((3/2 * cell.void * cell.u_B) - (1/2 * next_CFE.mesh[i-1].void * next_CFE.mesh[i-1].u_B) + ((-3/2 * next_CFE.mesh[i-1].void)+(1/2 * next_CFE.mesh[i-2].void * next_CFE.mesh[i-2].u_B))) + (cell.void * cell.cell_vol * cell.g * (cell.rho_H-cell.rho_U)) + cell.drag)
        
    Q_res = np.absolute(next_CFE.OG_MW*(10**6) - e_to_H2)
    mom_check = [mom,fld_mom,fld_arch,fld_drag,fld_cont,rho,fld_u]
    return [T_new,mom_check,next_CFE.X,Q_res,next_CFE.q_tot,rho_H, next_CFE.press_profile]

            


if __name__ =="__main__":

    PS_data_1 =np.array([-4871942.2439, 6862317.227, -3852036.5529, 1077768.9162, -150353.3467, 8370.8935])

    v_f = "1.1 - 23.333 * cell_radius" #"2.16 - 48 * cell_radius"
    T_P = "1500 "#+ 1.5*(20/(r))"
    q_tp = "10000000 + 150**(100*r)"
    BCs = ([1494,"dirichlet"],[0,"neumann"])

    base = {"P_core":10e6,
            "T_channel":450,
            "r5":56e-3,
            "d56":8e-3,
            "N":7000,
            "nu_s":0.691,
            "L_CFE":.84,
            "T_core":3700}

    stdlim = [0.5, 2]

    # ******************************************
    # TODO: if you don't want to plot vs a particular parameter, remove it from `vary`
    # TODO: if you want to plot a particular parameter over a range different from others, replace its limits in `vary`
    # ******************************************

    vary = {#"P_core":stdlim,
            #"T_channel":stdlim,
            # "r5":stdlim,
            # "d56":[0.125, 2],
            #"N":stdlim,
            #"nu_s":stdlim
            "L_CFE":stdlim
            }

    labels = {"P_core":"core pressure $P_3$", 
            "T_channel":"channel temperature $T_1$",
            "r5":"case radius $r_5$",
            "d56":"outer channel width $(r_6 - r_5)$",
            "N":"CFE rotation rate $\omega$",
            "nu_s":"turbine sp. speed $\\nu_s$",
            "L_CFE":"CFE length $l$"}

    yvals = dict()
    xvals = dict()

    base_core = H2(P=base["P_core"], T=base["T_core"])      # speed code up by not calculating on every single loop

    for key in vary:                    # iterate through properties we want to vary
        props = base.copy()             # reset all properties to base values
        lim = vary[key]                 # relative property value limits
        n_pts = 50                      # number of x-value points to plot to form a smooth curve
        xvals[key] = np.arange(lim[0], lim[1]+1e-5, np.diff(lim)/n_pts)
        yvals[key] = []
        i = 0                           # counter
        # pool = multi.Pool(15)
        CFEs = []
        results = []
        for x in xvals[key]:
            print(f"\n\tRELATIVE PROPERTY VALUE {key}: {x}")
            props[key] = x * base[key]              # adjust single parameter
            L_total = props["L_CFE"] + 0.1          # compute total length
            if key=="P_core" or key=="T_core":      # check if core state needs adjustment
                core = H2(P=props["P_core"], T=props["T_core"])
            else:
                core = base_core
            omega = props["N"] * np.pi/30           # compute omega
            r6 = props["r5"] + props["d56"]         # compute r6
            i +=1
            
            r1 = .03
            r2 = .045
            mdot = .108
            num_cells = 250
            P_core = props["P_core"]
            
            N = props["N"]
            L = props["L_CFE"]

            # CFEs.append()            
            # results = pool.map(solve, CFEs)
            # for index in range(len(results)):
            #R_1 = results[index]
            #D_1 = CFEs[index]
            
            D_1 = CFE(r2,r2-r1,L,7, num_cells, T_P, mdot, N, P_core,BCs,PS_data_1,OG_power=7)
            R_1 = solve(D_1)

            T_1 = R_1[0]
            M_1 = R_1[1]
            X_1 = R_1[2]
            QR_1 = R_1[3]
            Q_1 = R_1[4]
            
            x_1 = np.array(D_1.cell_center_radii)

            temp = T_1
            radius = np.array(x_1)
            void = np.array(X_1)
            fuel_density = np.array(M_1[5])
            prop_density = np.array(R_1[5])
            P_profile = R_1[6]

            mix_density = fuel_density * (1 - void) + prop_density * void
            np.save(f"Better/d_{key}_{i}.npy", mix_density)
            np.save(f"Better/r_{key}_{i}.npy", radius)
    # get baseline parameters
    base_cfe = CFE(r2,r2-r1,base["L_CFE"],7, num_cells, T_P, mdot, base["N"], base["P_core"],BCs,PS_data_1,OG_power=7)
    base_results = solve(base_cfe)
    T_1 = R_1[0]
    M_1 = R_1[1]
    X_1 = R_1[2]
    QR_1 = R_1[3]
    Q_1 = R_1[4]
    x_1 = np.array(base_cfe.cell_center_radii)
    temp = T_1
    radius = np.array(x_1)
    void = np.array(X_1)
    fuel_density = np.array(M_1[5])
    prop_density = np.array(R_1[5])
    P_profile = R_1[6]

    mix_density = fuel_density * (1 - void) + prop_density * void
    np.save(f"Better/d_baseline2.npy", mix_density)
    np.save(f"Better/r_baseline2.npy", radius)

    print("Design 1 Max Temp:",np.amax(T_1))
    print("Design 1 Max Void:",np.amax(X_1))
    print("Design 1 Q Deficit:",QR_1)
    print("Design 1 Scaled Q:",Q_1)

    plt.title("Density")
    plt.xlabel("Radius (m)")
    plt.ylabel("Density (kg/m^3)")
    plt.plot(x_1,M_1[5],label = "Uranium")
    plt.plot(x_1,prop_density,label = "Hydrogen")
    plt.plot(x_1,mix_density,label = "Mixture")
    plt.legend()
    plt.show()

    plt.title("Pressure")
    plt.xlabel("Radius (m)")
    plt.ylabel("Pressure (N/m^2)")
    plt.plot(x_1, P_profile, label = "Uranium")
    plt.legend()
    plt.show()
