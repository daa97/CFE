import numpy as np
np.set_printoptions(precision=1)
import matplotlib.pyplot as plt
import scipy as sp
import sys
from fluids import H2
import multiprocessing as multi
import prop_vary as pv


class Face:
    def __init__(self,face_radius, CFE_len):
        self.radius     = face_radius #m
        self.area       = 2 * np.pi * face_radius * CFE_len #m^2

class Cell:
    def __init__(self,cell_len,cell_temp,cell_radius,CFE_len,BCs,i,num_cells,press,m_b,mass_flow,omega): 
        self.press  = press
        self.T   = float(cell_temp) #FIXME - assign temp from temp profile
        self.H2     = H2(P=self.press, T=self.T, linear=True)
        self.rho_H  = self.H2.rho
        self.Cp     = self.H2.Cp
        self.h      = self.H2.h

        self.i           = i
        self.dR    = cell_len 
        self.BC          = BCs 
        self.Rm = cell_radius #Taken from CFE create_mesh method, mean radius

        self.omega  = omega #rpm
        self.m_dot  = mass_flow #kg/s

        self.k_U    = 13.7 #W/m-K 
        self.rho_U  = (19022 - 1.4486 * self.T) #kg/m^3
        self.Cp_U   = 179.78 + 0.0131 * self.T #J/kg/K - R**2 - 0.9998

        self.g      = self.omega ** 2 * self.Rm 
        self.V_b_o  = 1e-9 #m^3
        self.m_b    = m_b #kg
        self.m_b_o  = self.calc_init_bubble_mass() #kg
        self.V_B    = self.m_b_o * self.H2.volume #m^3  # individual bubble volume
        self.a_B    = ((3 * self.V_B )/ (4 * np.pi)) ** (1/3)       # equivalent spherical bubble radius
        self.r_B    = 2.179 * self.a_B                              # actual cap bubble radius
        self.u_B    = 2/3 * np.sqrt(self.r_B * self.g)              # average velocity  
        
        self.num_cells  = num_cells
        self.CFE_len    = CFE_len
        
        self.Ai = 2 * np.pi * CFE_len * (self.Rm-self.dR/2)         # interior face surface area
        self.Ao = 2 * np.pi * CFE_len * (self.Rm+self.dR/2)         # exterior face surface area


        self.vol   = np.pi * ((self.Rm+self.dR/2)**2 - (self.Rm-self.dR/2)**2) * self.CFE_len
        self.gas_vol    = self.m_dot * self.H2.volume * self.dR / self.u_B
        self.void       = self.calc_void_fraction()
        self.vol_void   = self.gas_vol / self.vol

        self.num_b      = self.calc_num_b()
        self.drag       = self.calc_drag()
        self.buoyant    = self.void * self.vol * self.g * (self.rho_H-self.rho_U)

    def calc_num_b(self):
        num_b = (self.m_dot/self.m_b_o) / self.u_B * self.dR
        return num_b

    def calc_void_fraction(self):
        N_B_dot = self.m_dot / self.Rm / self.u_B / self.rho_H
        void = N_B_dot / 2 / np.pi / self.CFE_len
        return void

    def calc_init_bubble_mass(self):
        if self.i ==0:
            m_b_o = self.rho_H * self.V_b_o
            return m_b_o
        else:
            return self.m_b

    def calc_drag(self):
        SA = np.pi * (self.r_B * np.sin(np.radians(52)))**2         # frontal area of moving bubble
        drag_per = 4/3 * self.rho_U * self.g * self.a_B * SA/self.V_B * self.void  
        drag = drag_per * self.vol
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
        self.dR = annulus_t/num_cells
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
        self.create_sys()

    def create_cell_center_radii(self):
        cell_center_radii = []
        for i in range(self.num_cells):
            cell_center_radii.append(self.OR_U-self.dR/2-self.dR*i)
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
            cell_vol = [(100**3)*cell.vol for cell in self.mesh]
            q_tp = np.multiply(q_t,cell_vol)
        return q_tp

    def find_press_profile(self):
        preint = (self.rpm * 2*np.pi/60)**2
        dPs = []; Ps = []
        rho_total = lambda cell: cell.rho_U * (1-cell.void) + cell.rho_H * cell.void
        for c in list(reversed(self.mesh)):
            dP = c.Rm * rho_total(c) * c.dR * preint
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
            mesh.append(Cell(cell_len = self.dR,cell_temp = self.temp_profile[i],cell_radius=self.cell_center_radii[i],CFE_len = self.L,
                        BCs = BCs,i=i,num_cells=self.num_cells,press=self.press_profile[i],m_b = m_b,mass_flow = self.mass_flow,omega = self.rpm*2 *np.pi/60))
            if i == 0:
                cell_0 = Cell(cell_len = self.dR,cell_temp = self.temp_profile[i],cell_radius=self.cell_center_radii[i],CFE_len = self.L,BCs = BCs,i=i,num_cells=self.num_cells,press=self.press_profile[i],m_b = m_b,mass_flow = self.mass_flow,omega = self.rpm*2 *np.pi/60)
                m_b = cell_0.m_b_o
        return mesh

    def create_sys(self): #Since Pe < 2 uses central differencing for discretization of convection term and a pseudo-transient for stability
        dt = 0.01
        q_tot = 0
        matrix = np.zeros((self.num_cells,self.num_cells),float)
        b = np.zeros((self.num_cells,1))
        X=[]
        for i,c in enumerate(self.mesh):
            X.append(c.void)

            a_c =(1-c.void) * c.k_U * (c.Ao + c.Ai) / c.dR + (1-c.void) * c.rho_U * c.Cp_U * c.vol / dt # diagonal matrix entry
            a_w = - (1-c.void) * c.k_U * c.Ao / c.dR                                                    
            a_e = - (1-c.void) * c.k_U * c.Ai / c.dR 

            SF_v = (1 - c.void)/.5
            q_tot += SF_v * self.q_tp[i]

            if i == 0 or i == (self.num_cells - 1):
                if c.BC[1].lower() == 'dirichlet':
                    if i == 0:
                        T_ghost =  2 * c.BC[0] - c.T #FIXME - Ghost Dirichlet BC is causing converging oscillations
                        ghost = Cell(self.dR,T_ghost,c.Rm+c.dR,self.L,BCs=[],i='',num_cells=self.num_cells,press=self.P_0,m_b = c.m_b_o,mass_flow = self.mass_flow,omega = self.rpm*2 *np.pi/60)
                        a_w_d =  - (1-c.void) * c.k_U * c.Ao / c.dR
                        matrix[i][i] = a_c - a_w_d
                        matrix[i][i+1] = a_e
                        b[i] = self.q_tp[i] * SF_v - 2 * a_w_d * c.BC[0] - self.mass_flow/2 * (self.mesh[i+1].h - ghost.h) + (1-c.void) * c.rho_U * c.Cp_U * c.vol * c.T / dt
                    else:
                        pass
                elif c.BC[1].lower() == 'neumann':
                    if i == 0:
                        pass
                    else:
                        a_e_n = - (1-c.void) * c.k_U * c.Ai / c.dR
                        matrix[i][i] = a_e_n + a_c
                        matrix[i][i-1] = a_w
                        b[i] = self.q_tp[i] * SF_v - a_e_n * c.BC[0] * c.dR - self.mass_flow/2 * (c.h - self.mesh[i-1].h) + (1-c.void) * c.rho_U * c.Cp_U * c.vol * c.T / dt
            else:
                matrix[i][i] = a_c
                matrix[i][i-1] = a_w
                matrix[i][i+1] = a_e
                b[i] = self.q_tp[i] * SF_v - self.mass_flow/2 * (self.mesh[i+1].h - self.mesh[i-1].h) + (1-c.void) * c.rho_U * c.Cp_U * c.vol * c.T / dt

        if np.abs(self.OG_MW*(10**6) - q_tot) >1e-4:
            mod = (self.OG_MW*(10**6))/q_tot
            self.MW = self.MW * mod
        
        self.matrix = matrix
        self.b = b
        self.q_tot = q_tot
        self.X = X

def status(msg):
    sys.stdout.write("\r"+msg+" "*20)
    sys.stdout.flush()

def solve(init_CFE):
    res = 100
    maxres = res
    tol = 1e-6
    c = 0
    old = init_CFE
    T_old = np.array(old.temp_profile).reshape(len(old.temp_profile),1)
    while  res > tol:
        if maxres == 100:           # get initial residual value
            maxres = res
        status(f"Err: {res/tol:.2f}; Completion: {(1-np.log(res/tol)/np.log(maxres/tol))*100:.3f}%")    # print progress 

        next = CFE(old.OR_U, old.annulus_t, old.L, old.MW, old.num_cells, T_old,old.mass_flow,
                        old.rpm, old.P_0,old.BC,old.q_profile_input,iteration = old.it + 1, OG_power = old.OG_MW, press_profile=old.find_press_profile())
        T_new = sp.linalg.solve(next.matrix,next.b)
        res = np.amax(np.divide(np.absolute(T_new-T_old),T_old))    

        T_old = T_new[:]
        old = next
        c += 1
        if c > 10000:
            res = 0
            
    next = CFE(old.OR_U,old.annulus_t, old.L, old.MW, old.num_cells, T_old,old.mass_flow, old.rpm, old.P_0, old.BC, old.q_profile_input,iteration = old.it + 1,OG_power=old.OG_MW, press_profile=old.find_press_profile())
    xi = np.array([next.BC[0][0], next.P_0])
    wall_h = H2(T=xi[0],P=xi[1], linear=True).h
    e_to_H2 = next.mass_flow * (next.mesh[next.num_cells-1].h - wall_h)

    T_ghost = 2 * next.mesh[0].BC[0] - next.mesh[0].T #FIXME - Ghost Dirichlet BC is causing converging oscillations
    ghost = Cell(next.dR,T_ghost,next.mesh[0].Rm+next.dR,next.L,BCs=[],i='',num_cells=next.num_cells,press=next.P_0,m_b = -1,mass_flow = next.mass_flow,omega = next.rpm*2 *np.pi/60)    
    mom=[]; fld_mom = []; fld_arch=[]; fld_drag=[]; fld_cont=[]; fld_u = []; rho = []; rho_H = []
    
    for i,c in enumerate(next.mesh):
        fld_u.append(c.u_B)
        rho.append(c.rho_U)
        rho_H.append(c.rho_H)
        fld_arch.append(c.buoyant)
        fld_drag.append(c.drag)

        if i == 0:
            cell_mom = (next.mass_flow / 2  * (next.mesh[i+1].void * next.mesh[i+1].u_B - ghost.void * ghost.u_B))
            fld_cont.append(1 / 2  * (next.mesh[i+1].void * next.mesh[i+1].rho_H * next.mesh[i+1].u_B * c.Ai   - ghost.void * ghost.rho_H * ghost.u_B * c.Ao))
            add = (next.mass_flow / 2  * (next.mesh[i+1].void * next.mesh[i+1].u_B - ghost.void * ghost.u_B)) + (c.void * c.vol * c.g * (c.rho_H-c.rho_U)) + c.drag
            mom.append(add)

        elif i == 1:
            cell_mom = next.mass_flow * ((3/2 * c.void * c.u_B) - (1/2 * next.mesh[i-1].void * next.mesh[i-1].u_B) + ((-3/2 * next.mesh[i-1].void) + (1/2 * ghost.void * ghost.u_B)))
            fld_cont.append(1/2 * (next.mesh[i+1].void * next.mesh[i+1].rho_H * next.mesh[i+1].u_B * c.Ai - next.mesh[i-1].void * next.mesh[i-1].rho_H * next.mesh[i-1].u_B*c.Ao))
            mom.append(next.mass_flow * ((3/2 * c.void * c.u_B) - (1/2 * next.mesh[i-1].void * next.mesh[i-1].u_B) + ((-3/2 * next.mesh[i-1].void) + (1/2 * ghost.void * ghost.u_B))) + (c.void * c.vol * c.g * (c.rho_H-c.rho_U))+ c.drag)
        
        elif i == next.num_cells-1:
            cell_mom = next.mass_flow * ((3/2 * c.void * c.u_B) - (1/2 * next.mesh[i-1].void * next.mesh[i-1].u_B) + ((-3/2 * next.mesh[i-1].void)+(1/2 * next.mesh[i-2].void * next.mesh[i-2].u_B))) 
            fld_cont.append(1/2 * (next.mesh[i].void * next.mesh[i].rho_H * next.mesh[i].u_B * c.Ai - next.mesh[i-1].void * next.mesh[i-1].rho_H * next.mesh[i-1].u_B*c.Ao))
            mom.append(next.mass_flow * ((3/2 * c.void * c.u_B) - (1/2 * next.mesh[i-1].void * next.mesh[i-1].u_B) + ((-3/2 * next.mesh[i-1].void)+(1/2 * next.mesh[i-2].void * next.mesh[i-2].u_B))) + (c.void * c.vol * c.g * (c.rho_H-c.rho_U)) + c.drag)
   
        else:
            cell_mom = next.mass_flow * ((3/2 * c.void * c.u_B) - (1/2 * next.mesh[i-1].void * next.mesh[i-1].u_B) + ((-3/2 * next.mesh[i-1].void)+(1/2 * next.mesh[i-2].void * next.mesh[i-2].u_B))) 
            fld_cont.append(1/2 * (next.mesh[i+1].void * next.mesh[i+1].rho_H * next.mesh[i+1].u_B * c.Ai - next.mesh[i-1].void * next.mesh[i-1].rho_H * next.mesh[i-1].u_B*c.Ao))
            mom.append(next.mass_flow * ((3/2 * c.void * c.u_B) - (1/2 * next.mesh[i-1].void * next.mesh[i-1].u_B) + ((-3/2 * next.mesh[i-1].void)+(1/2 * next.mesh[i-2].void * next.mesh[i-2].u_B))) + (c.void * c.vol * c.g * (c.rho_H-c.rho_U)) + c.drag)
        
        fld_mom.append(cell_mom)

    Q_res = np.absolute(next.OG_MW*(10**6) - e_to_H2)
    mom_check = [mom,fld_mom,fld_arch,fld_drag,fld_cont,rho,fld_u]
    return T_new,mom_check,next.X,Q_res,next.q_tot,rho_H, next.press_profile

def sweep(vary, num_cells,n_pts):        
    PS_data_1 =np.array([-4871942.2439, 6862317.227, -3852036.5529, 1077768.9162, -150353.3467, 8370.8935])
    v_f = "1.1 - 23.333 * cell_radius" #"2.16 - 48 * cell_radius"
    T_P = "1500 "#+ 1.5*(20/(r))"
    q_tp = "10000000 + 150**(100*r)"
    BCs = ([1494,"dirichlet"],[0,"neumann"])
    yvals = dict(); xvals = dict()
    
    base = pv.base
    for key in vary:                    # iterate through properties we want to vary
        props = base.copy()             # reset all properties to base values
        lim = vary[key]                 # relative property value limits
        xvals[key] = np.arange(lim[0], lim[1]+1e-5, np.diff(lim)/n_pts)
        yvals[key] = []
        i = 0                           # counter
        # pool = multi.Pool(15)
        results = []
        for x in xvals[key]:
            print(f"\n\tRELATIVE PROPERTY VALUE {key}: {x:.3f}")
            props[key] = x * base[key]              # adjust single parameter
            i +=1
            P_core = props["P_core"]
            N = props["N"]
            L = props["L_CFE"]
            r2 = props["r2"]
            r1 = props["r1"]

            cfe = CFE(r2,r2-r1,L,7, num_cells, T_P, props["mdot"], N, P_core,BCs,PS_data_1,OG_power=7)

            results = solve(cfe)
            temp = results[0]
            mom = results[1]
            void = np.array(results[2])
            radius = np.array(cfe.cell_center_radii)
            fuel_density = np.array(mom[5])
            prop_density = np.array(results[5])

            mix_density = fuel_density * (1 - void) + prop_density * void
            np.save(f"Better/d_{key}_{i}.npy", mix_density)
            np.save(f"Better/r_{key}_{i}.npy", radius)
            np.save(f"Better/v_{key}_{i}.npy", void)
            np.save(f"Better/t_{key}_{i}.npy", temp)
    
    # get baseline parameters
    N = base["N"]
    L = base["L_CFE"]
    r2 = base["r2"]
    r1 = base["r1"]
    base_cfe = CFE(r2,r2-r1,base["L_CFE"],7, num_cells, T_P, base["mdot"], base["N"], base["P_core"],BCs,PS_data_1,OG_power=7)
    results = solve(base_cfe)
    temp = results[0]
    mom = results[1]
    void = np.array(results[2])
    radius = np.array(base_cfe.cell_center_radii)
    
    fuel_density = np.array(mom[5])
    prop_density = np.array(results[5])

    mix_density = fuel_density * (1 - void) + prop_density * void
    np.save(f"Better/d_baseline.npy", mix_density)
    np.save(f"Better/r_baseline.npy", radius)
    np.save(f"Better/v_baseline.npy", void)
    np.save(f"Better/t_baseline.npy", temp)

    plt.title("Density")
    plt.xlabel("Radius (m)")
    plt.ylabel("Density (kg/m^3)")
    plt.plot(radius,fuel_density,label = "Uranium")
    plt.plot(radius,prop_density,label = "Hydrogen")
    plt.plot(radius,mix_density,label = "Mixture")
    plt.legend()
    plt.show()

if __name__ =="__main__":
    vary = pv.press_vary        # run all of the parameters

    
    
