import numpy as np
np.set_printoptions(precision=1)
from numpy.polynomial import Polynomial as poly
import matplotlib.pyplot as plt
import scipy as sp
import scipy.interpolate as inter
import csv
import pprint
import os
from timeit import Timer

#FIXME - Implement CUDA core solving of eqns 
# import numba
# import pyculib as pcl

#FIXME - General - replace large inputs with **kwargs

class Face:
    def __init__(self,face_radius, CFE_len):
        self.radius     = face_radius #m
        self.area       = 2 * np.pi * face_radius * CFE_len #m^2

class Cell:
    def __init__(self,cell_len,cell_temp,cell_radius,CFE_len,BCs,i,num_cells,press,m_b,mass_flow,omega): 
        self.i           = i
        self.cell_len    = cell_len 
        self.BC          = BCs 
        self.cell_radius = cell_radius #Taken from CFE create_mesh method
        self.cell_temp   = float(cell_temp) #FIXME - assign temp from temp profile
        self.R_H2        = 4124.2 #J/kg/K

        self.omega  = omega #rpm
        self.m_dot  = mass_flow #kg/s

        self.k_U    = 13.7 #W/m-K 
        self.rho_U  = self.calc_rho_U() #kg/m^3
        self.Cp_U   = 179.78 + 0.0131 * self.cell_temp #J/kg/K - R**2 - 0.9998
        self.press  = press*(10**6) #Pa

        self.V_b_o  = 1e-9 #m^3
        self.m_b    = m_b #kg
        self.m_b_o  = self.calc_init_bubble_mass() #kg

        self.V_B    = self.m_b_o * self.R_H2 * self.cell_temp / self.press #m^3
        self.u_B    = self.calc_avg_velocity() #m/s

        self.rho_H  = self.calc_rho_H2() #kg/m^3
        self.Cp     = self.get_Cp() #J/kg/K
        self.h      = self.get_h()  #J/kg

        self.num_cells  = num_cells
        self.CFE_len    = CFE_len
        self.g          = self.calc_g()
        
        self.face_e = Face(self.cell_radius-self.cell_len/2,self.CFE_len)
        self.face_w = Face(self.cell_radius+self.cell_len/2,self.CFE_len)

        self.cell_vol   = self.calc_cell_volume()
        self.gas_vol    = self.calc_gas_vol()
        self.void       = self.calc_void_fraction()
        self.vol_void   = self.calc_vol_void()
        self.num_b      = self.calc_num_b()
        self.drag       = self.calc_drag()

    def calc_avg_velocity(self):
        a_c = self.cell_radius * (self.omega ** 2)
        a = ((3 * self.V_B )/ (4 * np.pi)) ** (1/3)
        r_b = 2.179 * a
        u_B = 2/3 * np.sqrt(r_b * a_c)
        return u_B

    def calc_rho_H2(self):
        rho = self.press / (self.R_H2 * self.cell_temp)
        return rho
    
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
            rho_H = self.press / (4124.2 * self.cell_temp)
            m_b_o = rho_H * self.V_b_o
            return m_b_o
        else:
            return self.m_b

    def calc_cell_volume(self):
        cell_vol = np.pi * ((self.face_w.radius)**2 - (self.face_e.radius)**2) * self.CFE_len
        return cell_vol

    def calc_gas_vol(self):
        V_g = self.m_dot * self.R_H2 * self.cell_temp * self.cell_len / self.u_B / self.press
        return V_g

    def calc_vol_void(self):
        vol_void = self.gas_vol / self.cell_vol
        return vol_void

    def get_Cp(self):
        try:
            press = self.press / 101325 #kPa to atm
            points = (T_axis,P_axis)
            values = Cp_grid       
            xi = np.array([self.cell_temp, press])
            Cp = float(inter.interpn(points,values,xi)*1000)
            return Cp
        except ValueError:
            press = self.press / 101325 #kPa to atm
            points = (T_axis,P_axis)
            values = Cp_grid       
            xi = np.array([10000, press])
            Cp = float(inter.interpn(points,values,xi)*1000)
            return Cp

    def get_h(self):
        try:
            press = self.press / 101325 #kPa to atm
            points = (T_axis,P_axis)
            values = h_grid    
            xi = np.array([self.cell_temp, press])
            h = float(inter.interpn(points,values,xi)*1000)
            return h
        except ValueError:
            press = self.press / 101325 #kPa to atm
            points = (T_axis,P_axis)
            values = h_grid      
            if self.cell_temp > 0: 
                xi = np.array([10000, press])
                h = float(inter.interpn(points,values,xi)*1000)
            else:
                xi = np.array([400, press])
                h = float(inter.interpn(points,values,xi)*1000)
            return h

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
    def __init__(self,IR,annulus_t,length, power,num_cells,temp_profile,mass_flow,rpm,P_0,BC,q_tp,iteration=0,OG_power=0):
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
        self.mesh = self.create_mesh()
        self.q_profile_input = q_tp
        self.q_tp = self.create_heat_profile()
        self.sys = self.create_sys()
        self.matrix = self.sys[0]
        self.b = self.sys[1]
        self.q_tot = self.sys[2]
        self.X = self.sys[3]

        self.R_H2 = 4124.2 

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
            mesh.append(Cell(cell_len = self.cell_len,cell_temp = self.temp_profile[i],cell_radius=self.cell_center_radii[i],CFE_len = self.L,BCs = BCs,i=i,num_cells=self.num_cells,press=self.P_0,m_b = m_b,mass_flow = self.mass_flow,omega = self.rpm*2 *np.pi/60))
            if i == 0:
                cell_0 = Cell(cell_len = self.cell_len,cell_temp = self.temp_profile[i],cell_radius=self.cell_center_radii[i],CFE_len = self.L,BCs = BCs,i=i,num_cells=self.num_cells,press=self.P_0,m_b = m_b,mass_flow = self.mass_flow,omega = self.rpm*2 *np.pi/60)
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
                        notavar = ''
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
        else:
            pass
        print(sys)
        return [sys,b,q_tot,X]


def solve(init_CFE):
    res = 100
    tol = 1e-6
    c = 0
    old_CFE = init_CFE
    T_old = np.array(old_CFE.temp_profile).reshape(len(old_CFE.temp_profile),1)
    hist = np.zeros((init_CFE.num_cells-1,4))
    while  res > tol:

        next_CFE = CFE(old_CFE.OR_U,old_CFE.annulus_t,old_CFE.L, old_CFE.MW,old_CFE.num_cells,T_old,old_CFE.mass_flow,old_CFE.rpm,old_CFE.P_0,old_CFE.BC,old_CFE.q_profile_input,iteration = old_CFE.it + 1,OG_power = old_CFE.OG_MW)
        T_new = sp.linalg.solve(next_CFE.matrix,next_CFE.b)

        #FIXME - Implement CUDA-core computation with pyculib
        # desc = pcl.sparse.Sparse.matdescr()
        # R_nnz = np.linspace(next_CFE.num_cells)
        # nnz = pcl.sparse.Sparse.nnz('C',next_CFE.num_cells,next_CFE.num_cells,desc,next_CFE.matrix,R_nnz)
        # print(nnz)

        res = np.amax(np.divide(np.absolute(T_new-T_old),T_old))    

        T_old = T_new[:]
        old_CFE = next_CFE
        c += 1
        if c > 10000:
            res = 0          
            
    next_CFE = CFE(old_CFE.OR_U,old_CFE.annulus_t,old_CFE.L, old_CFE.MW,old_CFE.num_cells,T_old,old_CFE.mass_flow,old_CFE.rpm,old_CFE.P_0,old_CFE.BC,old_CFE.q_profile_input,iteration = old_CFE.it + 1,OG_power=old_CFE.OG_MW)
    points = (T_axis,P_axis)
    values = h_grid    
    xi = np.array([next_CFE.BC[0][0], next_CFE.P_0/0.101325])
    wall_h = float(inter.interpn(points,values,xi)*1000)
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
    for i,cell in enumerate(next_CFE.mesh):
        fld_u.append(cell.u_B)
        rho.append(cell.rho_U)
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

    print(e_to_H2)
    print(c)
    return [T_new,mom_check,next_CFE.X,Q_res,next_CFE.q_tot]

            
def create_interpolation_grid():
    Tmp = []
    Cp_grid =[]
    h_grid = []
    with open(os.path.join(os.path.dirname( __file__ ), 'H2_data_Cp.csv'),newline='') as csvfile:
        C_p = csv.reader(csvfile,csv.QUOTE_NONNUMERIC,delimiter=',')
        for i,row in enumerate(C_p):
            if i == 0:
                row.pop(0)
                P = [float(x) for x in row]

            else:
                Tmp.append(row[0])
                row.pop(0)
                Cp = [float(y) for y in row]
                Cp_grid.append(Cp)
    with open(os.path.join(os.path.dirname( __file__ ), 'H2_data_h.csv'),newline='') as csvfile:
        h_vals = csv.reader(csvfile,csv.QUOTE_NONNUMERIC,delimiter=',')

        for i,row in enumerate(h_vals):
            if i == 0:
                row.pop(0)
            else:
                row.pop(0)
                h = [float(k) for k in row]
                h_grid.append(h)
    T = [float(temp.replace(",","")) for temp in Tmp]


    return [T,P,Cp_grid,h_grid]

if __name__ =="__main__":

    PS2 = np.array([342.807, -6788.196, 53492.0865,-209274.7818, 406003.0282, -311150.776])
    PS3 = np.array([6023.7734,-137714.2642, 1257534.5563, -5731600.8451, 13036404.1890, -11834090.4699])

    PS_data_1 =np.array([-4871942.2439, 6862317.227, -3852036.5529, 1077768.9162, -150353.3467, 8370.8935])
    PS_data_2 = np.flip(PS2)
    PS_data_3 = np.flip(PS3)

    H2_grid_prop = create_interpolation_grid()

    Cp_grid = H2_grid_prop[2]
    T_axis = H2_grid_prop[0]
    P_axis = H2_grid_prop[1]
    h_grid = H2_grid_prop[3]
    v_f = "1.1 - 23.333 * cell_radius" #"2.16 - 48 * cell_radius"
    T_P = "1500 "#+ 1.5*(20/(r))"
    q_tp = "10000000 + 150**(100*r)"
    BCs = ([1494,"dirichlet"],[0,"neumann"])

    D_1 = CFE(4.5/100,(4.5-3)/100,84/100,7,10,T_P,0.108,7000,5,BCs,PS_data_1,OG_power=7)
    # D_2 = CFE(5.5/100,(5.5-3)/100,84/100,7,10,T_P,0.105,7000,10,BCs,PS_data_2,OG_power=7)
    # D_3 = CFE(5.5/100,(5.5-4)/100,84/100,11,10,T_P,0.155,7000,10,BCs,PS_data_3,OG_power=11)

    R_1 = solve(D_1)
    # R_2 = solve(D_2)
    # R_3 = solve(D_3)

    # T_1 = R_1[0]
    # T_2 = R_2[0]
    # T_3 = R_3[0]

    # M_1 = R_1[1]
    # M_2 = R_2[1]
    # M_3 = R_3[1]

    # X_1 = R_1[2]
    # X_2 = R_2[2]
    # X_3 = R_3[2]

    # QR_1 = R_1[3]
    # QR_2 = R_2[3]
    # QR_3 = R_3[3]  

    # Q_1 = R_1[4]
    # Q_2 = R_2[4]
    # Q_3 = R_3[4]  

    # x_1 = np.array(D_1.cell_center_radii)
    # x_2 = np.array(D_2.cell_center_radii)
    # x_3 = np.array(D_3.cell_center_radii)

    # print("Design 1 Max Temp:",np.amax(T_1))
    # print("Design 2 Max Temp:",np.amax(T_2))
    # print("Design 3 Max Temp:",np.amax(T_3))

    # print("Design 1 Max Void:",np.amax(X_1))
    # print("Design 2 Max Void:",np.amax(X_2))
    # print("Design 3 Max Void:",np.amax(X_3))

    # print("Design 1 Q Deficit:",QR_1)
    # print("Design 2 Q Deficit:",QR_2)
    # print("Design 3 Q Deficit:",QR_3)

    # print("Design 1 Scaled Q:",Q_1)
    # print("Design 2 Scaled Q:",Q_2)
    # print("Design 3 Scaled Q:",Q_3)

    # plt.title("Temperature Profile")
    # plt.xlabel("Radius (m)")
    # plt.ylabel("Temperature (K)")
    # plt.plot(x_1,T_1,label = "Design 1")
    # plt.plot(x_2,T_2,label = "Design 2")
    # plt.plot(x_3,T_3,label = "Design 3")
    # plt.legend()
    # plt.show()

    # plt.title("Void Fraction")
    # plt.xlabel("Radius (m)")
    # plt.ylabel("Void Fraction")
    # plt.plot(x_1,X_1,label = "Design 1")
    # plt.plot(x_2,X_2,label = "Design 2")
    # plt.plot(x_3,X_3,label = "Design 3")
    # plt.legend()
    # plt.show()

    # plt.title("Total Momentum Change Deficit")
    # plt.xlabel("Radius (m)")
    # plt.ylabel("Force (N)")
    # plt.plot(x_1,M_1[0],label = "Design 1")
    # plt.plot(x_2,M_2[0],label = "Design 2")
    # plt.plot(x_3,M_3[0],label = "Design 3")
    # plt.legend()
    # plt.show()

    # plt.title("Momentum Change Deficit - Linear Momentum Contribution")
    # plt.xlabel("Radius (m)")
    # plt.ylabel("Force (N)")
    # plt.plot(x_1[1:-1],M_1[1][1:-1],label = "Design 1")
    # plt.plot(x_2[1:-1],M_2[1][1:-1],label = "Design 2")
    # plt.plot(x_3[1:-1],M_3[1][1:-1],label = "Design 3")
    # plt.legend()
    # plt.show()

    # plt.title("Momentum Change Deficit - Archimedes Principle Contribution")
    # plt.xlabel("Radius (m)")
    # plt.ylabel("Force (N)")
    # plt.plot(x_1,M_1[2],label = "Design 1")
    # plt.plot(x_2,M_2[2],label = "Design 2")
    # plt.plot(x_3,M_3[2],label = "Design 3")
    # plt.legend()
    # plt.show()

    # plt.title("Momentum Change Deficit - Drag Contribution")
    # plt.xlabel("Radius (m)")
    # plt.ylabel("Force (N)")
    # plt.plot(x_1,M_1[3],label = "Design 1")
    # plt.plot(x_2,M_2[3],label = "Design 2")
    # plt.plot(x_3,M_3[3],label = "Design 3")
    # plt.legend()
    # plt.show()    

    # plt.title("Conservation of Mass")
    # plt.xlabel("Radius (m)")
    # plt.ylabel("Mass Flow Rate Imbalance")
    # plt.plot(x_1[0:-2],M_1[4][0:-2],label = "Design 1")
    # plt.plot(x_2[0:-2],M_2[4][0:-2],label = "Design 2")
    # plt.plot(x_3[0:-2],M_3[4][0:-2],label = "Design 3")
    # plt.legend()
    # plt.show() 

    # plt.title("Uranium Density")
    # plt.xlabel("Radius (m)")
    # plt.ylabel("Density (kg/m^3)")
    # plt.plot(x_1,M_1[5],label = "Design 1")
    # plt.plot(x_2,M_2[5],label = "Design 2")
    # plt.plot(x_3,M_3[5],label = "Design 3")
    # plt.legend()
    # plt.show() 

    # plt.title("Velocity Profile")
    # plt.xlabel("Radius (m)")
    # plt.ylabel("u (m/s)")
    # plt.plot(x_1,M_1[6],label = "Design 1")
    # plt.plot(x_2,M_2[6],label = "Design 2")
    # plt.plot(x_3,M_3[6],label = "Design 3")
    # plt.legend()
    # plt.show() 