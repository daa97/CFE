from pathlib import Path
import numpy as np
from icecream import ic
from matplotlib import pyplot as plt
import os
import subprocess as sp
import math
import sys

from Nozzledesigner.Rocket_CEA import create_CEA_infile
from Nozzledesigner.Ideal_rocket_perf import *
from fluids import *
import scipy.optimize
import scipy.integrate as integ

class Contour:
    def __init__(self, thrust=[],idl_gamma = [], idl_MW=[], p_c=[],p_amb=[],T_c=[],e=[],A_t=[],new_geometry=[],X=[],Y=[],showplot=False):
        #Must specify EITHER
        # 1. thrust, ideal specific heat ratio (idl_gamma), ideal molecular weight (idl_MW), chamber pressure (p_c),
        # ambient pressure (p_amb), chamber temperature (T_c), and exit expansion area ratio (e) OR

        # 2. X,Y (array of coordinates representing a predefined contour)
        # throat area (A_t), exit expansion area ratio (e)

        #optional parameters: new_geometry: change actual contour geometry (recommended for off-design). This needs optimization
        #showplot: self-explanatory


        # Option A: imported contour
        self.X = X
        self.Y = Y
        self.A_t = A_t

        default_geometry = {"L_star": .1, "theta_n": math.radians(30), "alpha": math.radians(15),
                            "theta_e": math.radians(8),
                            "cone_noz_perc": .8, "data_pts": 100, "e_c": 1.5}
        if new_geometry != []:
            for key, value in default_geometry.items():
                setattr(self, key, new_geometry.get(key, value))
        else:
            for key, value in default_geometry.items():
                setattr(self, key, default_geometry.get(key, value))


        # Needed for Option 1 and 2:
        self.e = e
        self.showplot = showplot

        # Option 2:
        if A_t==[]:
            self.A_t = self.ideal_properties(thrust=thrust, gamma=idl_gamma, p_c=p_c, p_amb=p_amb, T_c=T_c, e=e, MW=idl_MW)[0]
            self.make_contour()
        self.d_t = math.sqrt((self.A_t*4/math.pi))
        r_t=self.d_t/2
        self.r_c =1.5*r_t #Based on the specific approximation I am using to create the nozzle contour. This
        #approximation is listed both in H&H (Huzel and Huang) and in Sutton

    def ideal_properties(self, thrust,gamma,p_c,p_amb,T_c,e,MW):

        c_F_idl, c_star_idl, I_sp_idl, Tau_over_A_t, u_eq_idl, u_e_idl = rocket_perf(gamma, p_c, p_amb, T_c, 9.81, e, 8314.5/MW, [])
        A_t = thrust/Tau_over_A_t
        self.A_t = A_t
        self.e = e
        return A_t, c_F_idl, c_star_idl, I_sp_idl, u_eq_idl, u_e_idl

    def make_contour(self):

        #derived from inputs
        d_t =math.sqrt(self.A_t*4/math.pi)
        self.d_t = d_t
        r_t = d_t / 2
        r_e = np.sqrt(self.e)*r_t
        conv_cone_angle = math.radians(20)

        def cone_length(alpha, epsilon, r_t, r_e):
            return (r_t * (math.sqrt(self.e) - 1) + r_e * (1 / math.cos(self.alpha) - 1)) / math.tan(self.alpha)

        d_c = math.sqrt(self.e_c) * d_t
        r_c = d_c / 2
        conv_cone_len = cone_length(conv_cone_angle, self.e_c, r_t, r_c)
        L_n = self.cone_noz_perc * cone_length(conv_cone_angle, self.e, r_t, r_e)

        '''lmbda = (1+math.cos(self.alpha))/2

        def cone_length(alpha,epsilon,r_t,r_e):
            return (r_t * (math.sqrt(self.e) - 1) + r_e * (1/math.cos(self.alpha) - 1)) / math.tan(self.alpha)

        #Converging section
        V_c = self.L_star*self.A_t

        d_c = math.sqrt(self.e_c)*d_t

        r_c = d_c/2
        conv_cone_angle = math.radians(20)
        conv_cone_len = cone_length(conv_cone_angle, self.e_c,r_t,r_c)
        conv_cone_vol = math.pi/3*conv_cone_len*(r_c**2+r_t**2+r_c*r_t)

        V_c_cyl = V_c-conv_cone_vol
        L_c_cyl = V_c_cyl/(self.e_c*self.A_t)

        #Diverging section
        

        #N coordinates (point where it transitions to parabola
        N_t = .385*r_t*math.sin(math.radians(30))

d_c = math.sqrt(self.e_c)*d_t

        r_c = d_c/2
        N_a = r_t +.385*r_t*(1-math.cos(self.theta_n))
        E_t = L_n
        E_a = r_e'''


        # Region I

        r1 = 1.5*r_t
        th1 = np.linspace(-math.pi,-math.pi/2, int(math.floor(self.data_pts/3)))
        x1 = r1*np.cos(th1)
        y1 = r1*np.sin(th1)+2.5*r_t

        # Region II
        def n_point_finder(self,r_t):

            r = .385 * r_t
            slope = math.tan(self.theta_n)
            theta = np.linspace(-math.pi/2, 0, int(math.floor(self.data_pts/3)))
            theta.astype(int)
            x = r * np.cos(theta)
            y = r * np.sin(theta)+1.385*r_t
            grad_array = np.gradient(y, x)

            def find_nearest(array, value):
                idx = np.searchsorted(array, value, side="left")
                if idx > 0 and (idx == len(array) or math.fabs(value - array[idx - 1]) < math.fabs(value - array[idx])):
                    return idx - 1
                else:
                    return idx

            idx = find_nearest(grad_array, slope)
            x_n = x[idx]
            y_n = y[idx]
            return x_n, y_n


        x_n, y_n = n_point_finder(self, r_t)
        R = .385*r_t
        th2 = np.linspace(-np.pi/2,(self.theta_n-np.pi/2), int(math.floor(self.data_pts/3)))
        x2 = R*np.cos(th2)
        y2 = R*np.sin(th2)+1.385*r_t

        # Region III
        # fit parabola between two points
        y_e = r_e
        x_e = L_n
        m_n = math.tan(self.theta_n)
        m_e = math.tan(self.theta_e)

        A_mat = np.array([[2*y_n*m_n, m_n],
                     [2*y_e*m_e, m_e]])
        B_mat = np.array([1, 1])
        A,B = np.linalg.solve(A_mat, B_mat)
        C = x_n -A*y_n**2 -B*y_n

        y3 = np.array(np.linspace(y_n,y_e,int(math.floor(self.data_pts/3))))
        x3 = A*y3**2+B*y3+C

        # Plots

        if self.showplot:
            plt.plot(x1, y1, label='1')
            plt.plot(x2, y2, label='2')
            plt.plot(x3, y3, label='3')
            plt.axis('square')
            plt.show(block=False)
        X = np.concatenate((x1,x2,x3),axis=0)
        Y = np.concatenate((y1,y2,y3),axis=0)
        throat_ind = np.argmin(Y)

        subsonic_log= np.zeros((throat_ind),dtype=bool)

        supersonic_log = np.ones(((len(Y)-throat_ind)),dtype=bool)


        supersonic_check = np.concatenate((subsonic_log,supersonic_log),axis=0)
        self.X = X
        self.Y = Y
        self.supersonic_check = supersonic_check #This is an array of whether the flow is supersonic
        return X,Y, supersonic_check


class NozzleCEA:
    def __init__(self, contour, T_c,p_c,mdot):
        self.T_c = T_c
        self.p_c = p_c
        self.e = contour.e
        self.e_c = contour.e_c
        self.supersonic_check = contour.supersonic_check
        self.X = contour.X
        self.Y = contour.Y
        self.A_t = contour.A_t
        self.d_t = contour.d_t
        self.r_c = contour.r_c
        fzn_nz_props = {"P, BAR": [], "T, K": [], "RHO, KG/CU M": [], "H, KJ/KG": [], "U, KJ/KG": [],
                        "S, KJ/(KG)(K)": [],
                        "M, (1/n)": [], "GAMMAs": [], "SON VEL,M/SEC": [], "VISC,MILLIPOISE": [],
                        "MACH NUMBER": [], "Cp, KJ/(KG)(K)": [], "CONDUCTIVITY   ": [], "PRANDTL NUMBER": [],
                        "SON VEL,M/SEC": [], "*H":[], "*H2":[],"CSTAR, M/SEC":[]
                        }
        eql_nz_props = {"P, BAR": [], "T, K": [], "RHO, KG/CU M": [], "H, KJ/KG": [], "U, KJ/KG": [],
                        "S, KJ/(KG)(K)": [],
                        "M, (1/n)": [], "GAMMAs": [], "SON VEL,M/SEC": [], "VISC,MILLIPOISE": [],
                        "MACH NUMBER": [], "Cp, KJ/(KG)(K)": [], "CONDUCTIVITY   ": [], "PRANDTL NUMBER": [],
                        "SON VEL,M/SEC": [], "*H":[], "*H2":[],"CSTAR, M/SEC":[]
                        }
        prop_conv = {"P, BAR": 1e5, "T, K": 1, "RHO, KG/CU M": 1, "H, KJ/KG": 1e3, "U, KJ/KG": 1e3,
                     "S, KJ/(KG)(K)": 1e3, "M, (1/n)": 1, "GAMMAs": 1, "SON VEL,M/SEC": 1, "VISC,MILLIPOISE": 1e-4,
                     "MACH NUMBER": 1, "Cp, KJ/(KG)(K)": 1e3, "CONDUCTIVITY   ": 1, "PRANDTL NUMBER": 1,
                     "SON VEL,M/SEC": 1, "*H":1, "*H2":1, "CSTAR, M/SEC":1
                     }

        self.fzn_nz_props = fzn_nz_props
        self.eql_nz_props = eql_nz_props
        self.prop_conv = prop_conv
        areas = np.pi * np.array(self.Y) ** 2
        self.area_rats = areas / self.A_t
        if mdot:
            self.noz_heat = self.Noz_heat(mdot=mdot)


    def dict_updater(self, eql_nz_props,fzn_nz_props):
        for dictkey, dictval in eql_nz_props.items():
            if dictval:
                eqlvalue = dictval[0]*self.prop_conv[dictkey]

                self.eql_nz_props[dictkey].append(eqlvalue)
        for dictkey, dictval in fzn_nz_props.items():
            if dictval:
                fznvalue = dictval[0]* self.prop_conv[dictkey]
                self.fzn_nz_props[dictkey].append(fznvalue)

    def run_CEA(self,area_rat,supersonic):
        # Create CEA Input File
        from random import randint as rint
        os.chdir(r"C:\Users\tdham\OneDrive - Georgia Institute of Technology\Years\Third\Fall\NASA\CFE-Thermal-main")
        dir = "RCEAexec\\CEAexec-win\\"
        if supersonic == False:
            in_template = "rcea_sub.txt"
        else:
            in_template = "rcea_sup.txt"

        # Replace Input file  Area ratio, temperature, pressure, with relevant params
        with open(dir + in_template, mode='r') as f_in_template:
            contents = f_in_template.read()
            if supersonic == False:
                contents = contents.replace("<SUBSONIC_AREA_RATIO>", str(area_rat))
            else:
                contents = contents.replace("<SUPERSONIC_AREA_RATIO>", str(area_rat))
            contents = contents.replace("<CHAMBER_STAGNATION_TEMPERATURE>", str(self.T_c))
            contents = contents.replace("<CHAMBER_STAGNATION_PRESSURE>", str(self.p_c/6894.76))
            writename = f"Temp_RCEA_in{rint(0, 1000)}"

        with open(dir + writename + '.inp', mode='w') as f_in:
            f_in.write(contents)
        infile = writename
        os.chdir(r"C:\Users\tdham\OneDrive - Georgia Institute of Technology\Years\Third\Fall\NASA\CFE-Thermal-main\RCEAexec\CEAexec-win")

        # RUN CEA COMMAND:
        sp_status = sp.run(["echo", f"{infile}", "|", "FCEA2.exe"], stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
        cea_out_file = infile+".out"
        eql_nz_props, fzn_nz_props = self.CEA_output_parser(cea_out_file, option="a")

        # Update Nozzle Dictionary
        for dictkey, dictval in eql_nz_props.items():
            if dictval !=[]:
                eqlvalue = dictval[0]*self.prop_conv[dictkey]

                self.eql_nz_props[dictkey].append(eqlvalue)
        for dictkey, dictval in fzn_nz_props.items():
            if dictval != []:
                fznvalue = dictval[0]* self.prop_conv[dictkey]
                self.fzn_nz_props[dictkey].append(fznvalue)
        self.cea_in_file = infile+".inp"

    def CEA_output_parser(self,cea_out_file, option):

        if option == "c": #show combustion chamber props
            search_index = 0
        elif option == "t": #show throat properties
            search_index = 1
        elif option == "a": #show properties at a specified area ratio
            search_index = -1
        fzn_nz_props = {"P, BAR": [], "T, K": [], "RHO, KG/CU M": [], "H, KJ/KG": [], "U, KJ/KG": [],
                        "S, KJ/(KG)(K)": [],
                        "M, (1/n)": [], "GAMMAs": [], "SON VEL,M/SEC": [], "VISC,MILLIPOISE": [],
                        "MACH NUMBER": [], "Cp, KJ/(KG)(K)": [], "CONDUCTIVITY   ": [], "PRANDTL NUMBER": [],
                        "SON VEL,M/SEC": [], "*H":[], "*H2":[],"CSTAR, M/SEC":[]
                        }
        eql_nz_props = {"P, BAR": [], "T, K": [], "RHO, KG/CU M": [], "H, KJ/KG": [], "U, KJ/KG": [],
                        "S, KJ/(KG)(K)": [],
                        "M, (1/n)": [], "GAMMAs": [], "SON VEL,M/SEC": [], "VISC,MILLIPOISE": [],
                        "MACH NUMBER": [], "Cp, KJ/(KG)(K)": [], "CONDUCTIVITY   ": [], "PRANDTL NUMBER": [],
                        "SON VEL,M/SEC": [], "*H":[], "*H2":[],"CSTAR, M/SEC":[]
                        }
        equilibheader = "THEORETICAL ROCKET PERFORMANCE ASSUMING EQUILIBRIUM"
        frozenheader = "THEORETICAL ROCKET PERFORMANCE ASSUMING FROZEN COMPOSITION"
        equilib = False
        #ic(cea_out_file)
        with open(cea_out_file) as file:
            lines = file.readlines()

        # iterate over each line of the file
        for num in range(len(lines)):
            textline = lines[num]
            if equilibheader in textline:
                equilib = True
            elif frozenheader in textline:
                # print("stopped")
                equilib = False
            # go through properties
            for dictkey, dictvalue in eql_nz_props.items():
                # Check if the property is references in that line
                if dictkey in textline:
                    exit_prop =self.string_splitter(textline,search_index)
                    exit_prop_conv = exit_prop
                    if equilib:
                        eql_nz_props[dictkey].append(exit_prop_conv)
                    else:
                        fzn_nz_props[dictkey].append(exit_prop_conv)
            #ic(cea_out_file)
            self.cea_out_file = cea_out_file
        return eql_nz_props, fzn_nz_props

    def is_numeric(self,string):
        if '.' in string:
            string = string.replace('.', '')
        if '-' in string:
            string = string.replace('-', '')
        if string.isdecimal():
            numeric = True
        else:
            numeric = False
        return numeric

    def string_splitter(self,textline, search_index):
        split_string = textline.split()
        new_string = []
        for string in split_string:
            list_string = list(string)
            if '-' in string:
                for char in range(len(list_string)):
                    if list_string[char] == '-':
                        list_string.insert(char, ' ')
                        new = ''.join(str(e) for e in list_string)
                        new_string.append(new)
            else:
                new_string.append(string)
        split_string = ' '.join(str(e) for e in new_string)
        split_string = split_string.split()

        new_split_string = []
        for a in reversed(range(len(split_string))):

                if '-' in split_string[a] and self.is_numeric(split_string[a-1]):
                    new_value = float(split_string[a - 1]) * 10 ** float(split_string[a])
                    new_split_string.insert(0, new_value)

        if new_split_string != []:
            split_string = new_split_string
            exit_prop = split_string[search_index]
        else:
            exit_prop = float(split_string[search_index])
        return exit_prop

    def exec_Nozzle_CEA(self):
        for A in range(len(self.area_rats)):
            self.run_CEA(self.area_rats[A],self.supersonic_check[A])
            os.chdir(r"C:\Users\tdham\OneDrive - Georgia Institute of Technology\Years\Third\Fall\NASA\CFE-Thermal-main\RCEAexec\CEAexec-win")
            os.remove(self.cea_in_file)
            os.remove(self.cea_out_file)

        #Renaming to account for change in units (unit conversions were applied iteratively
        new_keys= ["P, PASCAL", "H, J/KG", "U, J/KG", "S, J/(KG)(K)", "VISC,N*SEC/ SQ M","Cp, J/(KG)(K)", "CONDUCTIVITY",]
        old_keys = ["P, BAR", "H, KJ/KG", "U, KJ/KG", "S, KJ/(KG)(K)", "VISC,MILLIPOISE","Cp, KJ/(KG)(K)", "CONDUCTIVITY   ",]
        def replace_key(newkey,oldkey,dictionary):
            d = dictionary
            d[newkey] = d[oldkey]
            del d[oldkey]
            return d

        fzn_nz_props2 = self.fzn_nz_props
        eql_nz_props2 = self.eql_nz_props
        for i in range(len(new_keys)):
            replace_key(new_keys[i],old_keys[i],fzn_nz_props2)
            replace_key(new_keys[i], old_keys[i], eql_nz_props2)

        #overwriting:
        self.eql_nz_props = eql_nz_props2
        self.fzn_nz_props = fzn_nz_props2

    def print_plots(self, frozen_check):
        if frozen_check:
            this_dict = self.fzn_nz_props
        else:
            this_dict = self.eql_nz_props
        for dictkey, dictval in this_dict.items():
            if dictval !=[]:

                plt.figure()
                plt.plot(self.X,dictval,'bo',label=dictkey)
                plt.draw()
                plt.title(dictkey)
                plt.show(block=False)

    class Noz_heat:
        def __init__(self,mdot):
            w = .0254*1/8 #channel width (m)
            h = .0254*1/8 #channel depth (m)
            g = .0254*1/8 # maximum gap between channels at throat area (m)
            self.default_chan_geo = {"k_w":398,"t_w":.000254,"chan_width": w,"chan_depth":h, "chan_gap":g}
            self.mdot=mdot

            # using 1/8" square channel as default but NEEDS OPTIMIZATION
            # Using wall thickness of 10 thou as default but NEEDS OPTIMIZATION
            # Using copper as nozzle material as default but NEEDS OPTIMIZATION

        def compute_heatXfer(self, Nozzle, dist, nz_coolant_start,showplots,chan_geo):
            #Unpacking hot gas property Distributions

            T_hg_dist = np.array(dist["T, K"])
            rho_hg_dist = np.array(dist["RHO, KG/CU M"])
            sonic_vel_dist = np.array(dist["SON VEL,M/SEC"])
            M_dist = np.array(dist["MACH NUMBER"])
            vel_hg_dist = sonic_vel_dist*M_dist
            mu_hg_dist = np.array(dist["VISC,N*SEC/ SQ M"])
            cp_hg_dist = np.array(dist["Cp, J/(KG)(K)"])
            k_hg_dist = np.array(dist["CONDUCTIVITY"])
            Pr_hg_dist = np.array(dist["PRANDTL NUMBER"])
            cstar_dist = np.array(dist["CSTAR, M/SEC"])
            gamma_dist = np.array(dist["GAMMAs"])

            if chan_geo ==[]:
                chan_geo = self.default_chan_geo

            #Cooling Channel Geometry:
            k_w= chan_geo["k_w"] #W/m*K
            t_w = chan_geo["t_w"]
            char_len = t_w
            h = chan_geo["chan_depth"]
            w = chan_geo["chan_width"]
            g = chan_geo["chan_gap"]
            A_c = h*w       #Channel area
            N_c = math.floor(math.pi*Nozzle.d_t/(w+g)) #Number of channels

            # Convection Coefficients Estimation:
            # hot-gas side convective coefficient: h_g
            # coolant-side convective coefficient: h_L (coolant will start as liquid then become supercritical gas)

            # Coolant-Side Heat Transfer Coefficients

            '''nz_coolant = nz_coolant_start  # assuming one state
            k_L = nz_coolant.k
            rho_L = nz_coolant.rho
            v_L = self.mdot / (rho_L * A_c)
            mu_L = nz_coolant.mu
            Re = rho_L * v_L * char_len / mu_L
            cp_L = nz_coolant.cp
            Pr = mu_L * cp_L / k_L
            h_L = k_L / char_len * .023 * Re ** .8 * Pr ** .4  # Dittus Boelter correlation for turblent flow pipes
            T_L = nz_coolant.t  # assume constant for now
            h_L1 = h_L
            T_L1 = T_L

            # Hot-gas side heat transfer coefficients:
            # Low-Fidelity Approximation from H&H (It uses imperial units)

            rho_hg_dist_imperial = rho_hg_dist*2.20462/61023.7 #Converts from kg/m^3 to lb/in^3
            vel_hg_dist_imperial = vel_hg_dist *39.3701 #m/s to in/s
            h_g_dist_lf_imperial = (rho_hg_dist_imperial * vel_hg_dist_imperial) ** .8
            h_g_dist_lf2 = h_g_dist_lf_imperial* 1.635 * 10 ** 6 *1.8 #from BTU/(in^2*s*F) to W/(m^2*K)

            h_g_dist_lf = (rho_hg_dist * vel_hg_dist) ** .8
            # Even though this is not what H&H reads (H&H says you should use imperial), this produces an estimate closer to Bartz, so I am taking this to represent the low-fidelity
            #calculation for h_g
            # Heat Transfer Using Low-Fidelity Coefficients:
            q_dot_r = 0
            q_dist_lf = (T_hg_dist-T_L +q_dot_r/h_g_dist_lf)/(1/h_g_dist_lf +t_w/k_w + 1/h_L)
            q_dot_total_lf = N_c*integ.trapezoid(y=q_dist_lf,x=Nozzle.X)*w #multiply by w to get in units of W from W/m
            #( this assumes that the heat transfer is 1D and only is going through the channel area that is touching the
            # chamber; the integration across the nozzle contour is one length dimension, and the channel dimension is the other
            #N_c is the number of channels

            self.q_dot_total_lf = q_dot_total_lf'''


            # Higher-fidelity method: Bartz correlation:
            # This method is copied from the following technical note by Bartz in 1957:"A Simple Equation for Rapid Estimation of Rocket Nozzle Convective
            # Heat Transfer Coefficients "
            # http://arc.aiaa.org | DOI: 10.2514/8.12572

            def bartz(d_t, mu_c, Cp_c, Pr_c, p_c, T_c, Cstar, r_c, A, gamma, M, T_w,
                      w=.6):  # SI METRIC INPUTS (N,Pa,kg,s,K)

                # note: mu and Cp are considered constant, but I may be able to compute it at every step later on

                # r_c is the radius of curvature of the throat
                # A is the area where the heat transfer coefficient needs to be calculated
                # w is the correlation of mu with Temperature (mu ~ T^w) From Huzel and Huang, it is assumed to be .6

                # Converting metric to Imperial (equation is in imperial)
                A_t = math.pi / 4 * d_t ** 2  # should be in m^2 because d_t comes in as meters
                # A_t doesn't need to be converted because A and A_t are divided, and they are both the same unit
                d_t = d_t * 39.3701  # meters to inches
                mu_c = mu_c * 2.20462 / 39.3701  # kg/(m*s) to lbm/(in*s)
                Cp_c = Cp_c * 0.000947817 / 2.20462 / 1.8  # J/kgK to BTU/(lbm*F)
                p_c = p_c / 6894.76  # Pa to psi
                T_c = T_c * 9 / 5  # Kelvin to Rankine
                Cstar = Cstar * 3.28084  # m/s to fps
                r_c = r_c * 39.3701  # m to inches
                T_w = T_w * 9 / 5  # Kelvin to Rankine

                sigma = ((.5 * T_w / T_c * (1 + (gamma - 1) / 2 * M ** 2) + .5) ** (.8 - (w / 5)) * (
                        1 + (gamma - 1) / 2 * M ** 2) ** (w / 5)) ** -1
                h_g = (0.026 / (d_t ** .2) * (mu_c ** .2 * Cp_c / (Pr_c ** .6)) * (p_c * 32.2 / Cstar) ** .8 * (
                        d_t / r_c) ** .1) * (A_t / A) ** .9 * sigma

                h_g = h_g * 1.635 * 10 ** 6 * 1.8  # BTU/(in^2 sec F) to N/(m^2 sec K)
                return h_g  # Ns/m^2

            h_g_dist_hf = []
            q_dist_hf = []
            T_L_dist = []
            h_L_dist = []
            T_wg_dist = []

            for i in range(len(Nozzle.area_rats)):
                # Coolant-Side Heat Transfer Coefficients (needed to solve iteratively)
                T_hg = T_hg_dist[i]
                mu_hg = mu_hg_dist[i]
                cp_hg = cp_hg_dist[i]
                Pr_hg= Pr_hg_dist[i]
                gamma = gamma_dist[i]
                M=M_dist[i]
                cstar =cstar_dist[i]
                if i ==0:
                    nz_coolant = nz_coolant_start  # assuming one state
                else:
                    if not ("Updated" in os.getcwd()):
                        os.chdir(
                            r"C:\Users\tdham\OneDrive - Georgia Institute of Technology\Years\Third\Fall\NASA\CFE-Thermal-main")
                    H2 = Fluid("Hydrogen", prop_files)
                    nz_coolant= H2(T=q_dot/(self.mdot*cp_L)+T_L, P=700 * 6894.76)   # Assume no pressure drop in channels for now # CHECK IF THIS ENERGY EQUATION IS VALID IS VALID

                k_L = nz_coolant.k
                rho_L = nz_coolant.rho
                v_L = self.mdot / (rho_L * A_c)
                mu_L = nz_coolant.mu
                Re = rho_L * v_L * char_len / mu_L
                cp_L = nz_coolant.cp
                Pr = mu_L * cp_L / k_L

                h_L = k_L / char_len * .023 * Re ** .8 * Pr ** .4  # Dittus Boelter correlation for turblent flow pipes
                T_L = nz_coolant.t

                T_L_dist = np.append(T_L_dist,T_L)
                h_L_dist = np.append(h_L_dist,h_L)

                def eq2iteratethrough(T_wg):
                    ## GO BACK AND FIND OUT WHAT w should BE!!
                    # This equation below comes from assuming 1-D steady-state heat transfer with negligible radiation heat transfer
                    h_g = bartz(d_t=Nozzle.d_t, p_c =Nozzle.p_c,T_c = Nozzle.T_c, r_c=Nozzle.r_c,A=Nozzle.area_rats[i], mu_c=mu_hg, Cp_c =cp_hg,Pr_c = Pr_hg,
                          Cstar=cstar, gamma=gamma, M=M, T_w=T_wg)
                    eq= -T_wg +T_hg - (T_hg-T_L)/(1+(t_w/k_w +1/h_L)*h_g)

                    return eq

                T_wg = scipy.optimize.brentq(eq2iteratethrough, 0,Nozzle.T_c+100) #Solves nonlinear equation
                T_wg_dist = np.append(T_wg_dist, T_wg)

                h_g = bartz(d_t=Nozzle.d_t, p_c=Nozzle.p_c, T_c=Nozzle.T_c, r_c=Nozzle.r_c, A=Nozzle.area_rats[i],
                            mu_c=mu_hg, Cp_c=cp_hg, Pr_c=Pr_hg,
                            Cstar=cstar, gamma=gamma, M=M, T_w=T_wg)

                q_dot_r = 0
                q_dot = (T_hg - T_L + q_dot_r / h_g) / (1 / h_g + t_w / k_w + 1 / h_L)
                q_dist_hf = np.append(q_dist_hf,q_dot)
                h_g_dist_hf = np.append(h_g_dist_hf, h_g)


            # Heat TransferUsing High-Fidelity Coefficients:

            q_dot_total_hf = N_c*integ.trapezoid(y=q_dist_hf,x=Nozzle.X)*w  #multiply by w to get in units of W from W/m
            #( this assumes that the heat transfer is 1D and only is going through the channel area that is touching the
            # chamber; the integration across the nozzle contour is one length dimension, and the channel dimension is the other

            self.q_dot_total_hf = q_dot_total_hf

            # Plots
            if showplots:
                title = "Nozzle Heat Transfer/Temperature Distribution"
                fig, (ax1,ax2) = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=False)
                ax1.plot(Nozzle.X, T_hg_dist, 'bo',label="Temp (K)")
                ax1.set(xlabel="Axial Distance (X) (m)")
                ax1.set(ylabel="Temperature (K)")

                ax2.set(ylabel="Heat Flux Rate (W/m^2)")
                ax2.plot(Nozzle.X, q_dist_hf, 'ko',label="Heat Flux Rate (W/m^2) (high-fidelity Bartz)")
                #ax2.plot(Nozzle.X, q_dist_lf, 'ro', label="Heat Flux Rate (W/m^2) (Low-fidelity)")
                fig.suptitle(title)
                plt.grid(which='both')
                plt.show(block=False)
                plt.legend()
                ic(q_dot_total_hf)
                #ic(q_dot_total_lf)
            return q_dot_total_hf

        def CEA(self,Nozzle,frozen,nz_coolant_start,showplots,chan_geo=[]): # takes in properties dictionary and outputs heat transfer
            if frozen:
                property_distributions = Nozzle.fzn_nz_props
            else:
                property_distributions = Nozzle.eql_nz_props
            q_regen = self.compute_heatXfer(Nozzle, property_distributions, nz_coolant_start, showplots, chan_geo)
            return q_regen

        def oneDIsentropic(self,Nozzle,showplots):
            MW = 2 # for now
            R = 8314/MW #J/kgK

            area_rats = Nozzle.area_rats
            gamma_c = Nozzle.eql_nz_props["GAMMAs"][0]
            rho_c = Nozzle.eql_nz_props["RHO, KG/CU M"][0]

            if not ("Updated" in os.getcwd()):
                os.chdir(r"C:\Users\tdham\OneDrive - Georgia Institute of Technology\Years\Third\Fall\NASA\CFE-Thermal-main")
            H2 = Fluid("Hydrogen", prop_files)

            chamber_props = H2(T=Nozzle.T_c, P=Nozzle.p_c)
            gamma2 = chamber_props.Y
            ic(gamma_c)
            ic(gamma2)
            rho2 = chamber_props.rho
            ic(rho_c)
            ic(rho2)
            # Create array of mach numbers from area ratios
            M_subsonic = np.array([])
            M_supersonic = np.array([])
            throat_ind = np.argmin(area_rats)
            for i in area_rats[0:throat_ind]:
                M_sub = AonAStar.mach(self, Arat=i, gamma=gamma_c)[0]
                if M_sub.size ==0:
                    M_sub = M_subsonic[-1]
                M_subsonic = np.append(M_subsonic, M_sub)

            for i in area_rats[throat_ind:(len(area_rats))]:
                M_super = AonAStar.mach(self, Arat=i, gamma=gamma_c)[1]
                if M_super.size ==0:
                    M_super= M_supersonic[-1]
                M_supersonic = np.append(M_supersonic, M_super)


            M_dist = np.concatenate((M_subsonic, M_supersonic))

            plt.plot(Nozzle.area_rats,M_dist)
            plt.show(block=False)
            T_dist = []
            V_dist = []
            rho_dist = []
            M_dist = np.array(M_dist)

            # create 1D thermo properties from mach number distribution
            for i in M_dist:
                TO_over_T = compflowtool(gamma_c, i)[2]
                rho_O_over_rho = compflowtool(gamma_c, i)[3]

                T = Nozzle.T_c / TO_over_T
                V = np.sqrt(gamma_c * R * T) * i
                rho = rho_c / rho_O_over_rho

                T_dist = np.append(T_dist, T)
                V_dist = np.append(V_dist, V)
                rho_dist = np.append(rho_dist,rho)

            T_dist = np.array(T_dist)
            V_dist = np.array(V_dist)
            rho_dist = np.array(rho_dist)


            T_dist_CEA_eql = np.array(Nozzle.eql_nz_props["T, K"])
            rho_dist_CEA_eql = np.array(Nozzle.eql_nz_props["RHO, KG/CU M"])
            V_dist_CEA_eql = np.array(Nozzle.eql_nz_props["SON VEL,M/SEC"])*np.array(Nozzle.eql_nz_props["MACH NUMBER"])
            M_dist_CEA_eql = np.array(Nozzle.eql_nz_props["MACH NUMBER"])

            T_dist_CEA_fzn = np.array(Nozzle.fzn_nz_props["T, K"])
            rho_dist_CEA_fzn = np.array(Nozzle.fzn_nz_props["RHO, KG/CU M"])
            V_dist_CEA_fzn = np.array(Nozzle.fzn_nz_props["SON VEL,M/SEC"]) * np.array(Nozzle.fzn_nz_props["MACH NUMBER"])
            M_dist_CEA_fzn = np.array(Nozzle.fzn_nz_props["MACH NUMBER"])
            if showplots:
                fig, axs = plt.subplots(4)
                fig.suptitle('1-D Isentropic Flow versus NASA CEA 1-D Isentropic Flow')

                axs[0].set(ylabel="Mach Number")
                axs[0].plot(Nozzle.X, M_dist, '-ko', label="1-D")
                axs[0].plot(Nozzle.X, M_dist_CEA_eql, '--r*', label="CEA, Shifting")
                # axs[0].plot(Nozzle.X, M_dist_CEA_fzn, '--b*', label="CEA, Frozen")
                plt.grid(which='both')

                axs[1].set(ylabel="Temperature (K)")
                axs[1].plot(Nozzle.X, T_dist, '-ko',label="1-D")
                axs[1].plot(Nozzle.X, T_dist_CEA_eql, '--r*', label="CEA, Shifting")
                #axs[1].plot(Nozzle.X, T_dist_CEA_fzn, '--b', label="CEA, Frozen")
                plt.grid(which='both')

                axs[2].set(ylabel="Velocity (m/s)")
                axs[2].plot(Nozzle.X, V_dist, '-ko',label="1-D")
                axs[2].plot(Nozzle.X, V_dist_CEA_eql, '--r*', label="CEA, Shifting")
                #axs[2].plot(Nozzle.X, V_dist_CEA_fzn, '--b*', label="CEA, Frozen")
                plt.grid(which='both')

                axs[3].set(ylabel="Density (kg/m^3)")
                axs[3].plot(Nozzle.X, rho_dist,'-ko',label="1-D")
                axs[3].plot(Nozzle.X, rho_dist_CEA_eql, '--r*', label="CEA, Shifting")
                #axs[3].plot(Nozzle.X, rho_dist_CEA_fzn, '--b*', label="CEA, Frozen")
                plt.grid(which='both')
                plt.legend()

                plt.show(block=False)


# =====================================THERMODYNAMIC CYCLE CODE=================================================
if not ("Updated" in os.getcwd()):
    os.chdir(r"C:\Users\tdham\OneDrive - Georgia Institute of Technology\Years\Third\Fall\NASA\CFE-Thermal-main")
H2 = Fluid("Hydrogen", prop_files)

default_w = .0254 * 1 / 8  # channel width (m)
default_h = .0254 * 1 / 8  # channel depth (m)
default_g = .0254 * 1 / 8  # gap between channels (m)

default_chan_geo = {"k_w": 398, "t_w": .00254, "chan_width": default_w, "chan_depth": default_h, "chan_gap":default_g}

def get_q_nozzle(p_c=550*6894.76,T_c=3800,thrust = 10e3,e=120, mdot = 2.3, chan_geo=default_chan_geo, nz_cool_start = H2(T=120,P=700*6894.76), data_pts=25, frozen=False):
    chamber=H2(P=p_c,T=T_c)
    this_Contour = Contour(thrust=thrust,p_c=p_c,p_amb=0,idl_gamma=chamber.y,idl_MW=2, T_c=T_c,e=e, showplot=False, new_geometry={"data_pts": data_pts})
    Nozzle = NozzleCEA(this_Contour,T_c=T_c,p_c=p_c,mdot=mdot)
    Nozzle.exec_Nozzle_CEA()
    q_regen = Nozzle.noz_heat.CEA(Nozzle, nz_coolant_start=nz_cool_start, frozen=frozen, showplots=False, chan_geo=chan_geo)
    return q_regen




