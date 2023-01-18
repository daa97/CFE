#This code will allow you to input params into NASA CEA and run it for a given nozzle, given area ratio (no contour) frozen and equilibrium

from icecream import ic
from matplotlib import pyplot as plt
import os
import subprocess as sp
from fluids import *

class NozzleCEA:
    def __init__(self, T_c, p_c):
        self.T_c = T_c
        self.p_c = p_c

        fzn_nz_props = {"P, BAR": [], "T, K": [], "RHO, KG/CU M": [], "H, KJ/KG": [], "U, KJ/KG": [],
                        "S, KJ/(KG)(K)": [],
                        "M, (1/n)": [], "GAMMAs": [], "SON VEL,M/SEC": [], "VISC,MILLIPOISE": [],
                        "MACH NUMBER": [], "Cp, KJ/(KG)(K)": [], "CONDUCTIVITY   ": [], "PRANDTL NUMBER": [],
                        "SON VEL,M/SEC": [], "*H": [], "*H2": []
                        }
        eql_nz_props = {"P, BAR": [], "T, K": [], "RHO, KG/CU M": [], "H, KJ/KG": [], "U, KJ/KG": [],
                        "S, KJ/(KG)(K)": [],
                        "M, (1/n)": [], "GAMMAs": [], "SON VEL,M/SEC": [], "VISC,MILLIPOISE": [],
                        "MACH NUMBER": [], "Cp, KJ/(KG)(K)": [], "CONDUCTIVITY   ": [], "PRANDTL NUMBER": [],
                        "SON VEL,M/SEC": [], "*H": [], "*H2": []
                        }
        prop_conv = {"P, BAR": 1e5, "T, K": 1, "RHO, KG/CU M": 1, "H, KJ/KG": 1e3, "U, KJ/KG": 1e3,
                     "S, KJ/(KG)(K)": 1e3, "M, (1/n)": 1, "GAMMAs": 1, "SON VEL,M/SEC": 1, "VISC,MILLIPOISE": 1e-4,
                     "MACH NUMBER": 1, "Cp, KJ/(KG)(K)": 1e3, "CONDUCTIVITY   ": 1, "PRANDTL NUMBER": 1,
                     "SON VEL,M/SEC": 1, "*H": 1, "*H2": 1
                     }
        self.fzn_nz_props = fzn_nz_props
        self.eql_nz_props = eql_nz_props
        self.prop_conv = prop_conv
        areas = np.pi * np.array(self.Y) ** 2
        self.area_rats = areas / self.A_t
        self.noz_heat = self.Noz_heat()

    def run_CEA(self, area_rat, supersonic):
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
            contents = contents.replace("<CHAMBER_STAGNATION_PRESSURE>", str(self.p_c))
            writename = f"Temp_RCEA_in{rint(0, 1000)}"

        with open(dir + writename + '.inp', mode='w') as f_in:
            f_in.write(contents)
        infile = writename
        os.chdir(r"C:\Users\tdham\OneDrive - Georgia Institute of Technology\Years\Third\Fall\NASA\CFE-Thermal-main\RCEAexec\CEAexec-win")

        # RUN CEA COMMAND:
        sp_status = sp.run(["echo", f"{infile}", "|", "FCEA2.exe"], stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
        cea_out_file = infile + ".out"
        eql_nz_props, fzn_nz_props = self.CEA_output_parser(cea_out_file, option="a")

        # Update Nozzle Dictionary
        for dictkey, dictval in eql_nz_props.items():
            if dictval != []:
                eqlvalue = dictval[0] * self.prop_conv[dictkey]

                self.eql_nz_props[dictkey].append(eqlvalue)
        for dictkey, dictval in fzn_nz_props.items():
            if dictval != []:
                fznvalue = dictval[0] * self.prop_conv[dictkey]
                self.fzn_nz_props[dictkey].append(fznvalue)

    def CEA_output_parser(self, cea_out_file, option):

        if option == "c":  # show combustion chamber props
            search_index = 0
        elif option == "t":  # show throat properties
            search_index = 1
        elif option == "a":  # show properties at a specified area ratio
            search_index = -1
        fzn_nz_props = {"P, BAR": [], "T, K": [], "RHO, KG/CU M": [], "H, KJ/KG": [], "U, KJ/KG": [],
                        "S, KJ/(KG)(K)": [],
                        "M, (1/n)": [], "GAMMAs": [], "SON VEL,M/SEC": [], "VISC,MILLIPOISE": [],
                        "MACH NUMBER": [], "Cp, KJ/(KG)(K)": [], "CONDUCTIVITY   ": [], "PRANDTL NUMBER": []
                        }
        eql_nz_props = {"P, BAR": [], "T, K": [], "RHO, KG/CU M": [], "H, KJ/KG": [], "U, KJ/KG": [],
                        "S, KJ/(KG)(K)": [], "M, (1/n)": [], "GAMMAs": [], "SON VEL,M/SEC": [], "VISC,MILLIPOISE": [],
                        "MACH NUMBER": [], "Cp, KJ/(KG)(K)": [], "CONDUCTIVITY   ": [], "PRANDTL NUMBER": []
                        }
        equilibheader = "THEORETICAL ROCKET PERFORMANCE ASSUMING EQUILIBRIUM"
        frozenheader = "THEORETICAL ROCKET PERFORMANCE ASSUMING FROZEN COMPOSITION"
        equilib = False
        # ic(cea_out_file)
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
                    exit_prop = self.string_splitter(textline, search_index)
                    exit_prop_conv = exit_prop
                    if equilib:
                        eql_nz_props[dictkey].append(exit_prop_conv)
                    else:
                        fzn_nz_props[dictkey].append(exit_prop_conv)
            # ic(cea_out_file)
            self.cea_out_file = cea_out_file
        return eql_nz_props, fzn_nz_props

    def is_numeric(self, string):
        if '.' in string:
            string = string.replace('.', '')
        if '-' in string:
            string = string.replace('-', '')
        if string.isdecimal():
            numeric = True
        else:
            numeric = False
        return numeric

    def string_splitter(self, textline, search_index):
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

            if '-' in split_string[a] and self.is_numeric(split_string[a - 1]):
                new_value = float(split_string[a - 1]) * 10 ** float(split_string[a])
                new_split_string.insert(0, new_value)

        if new_split_string != []:
            split_string = new_split_string
            exit_prop = split_string[search_index]
        else:
            exit_prop = float(split_string[search_index])
        return exit_prop

    def print_plots(self, frozen_check):
        if frozen_check:
            this_dict = self.fzn_nz_props
        else:
            this_dict = self.eql_nz_props
        for dictkey, dictval in this_dict.items():
            if dictval != []:
                plt.figure()
                plt.plot(self.X, dictval, 'bo', label=dictkey)
                plt.draw()
                plt.title(dictkey)
                plt.show(block=False)

# Inputs:
if not ("Updated" in os.getcwd()):
    os.chdir(r"C:\Users\tdham\OneDrive - Georgia Institute of Technology\Years\Third\Fall\NASA\CFE-Thermal-main")
H2 = Fluid("Hydrogen", prop_files)
p_c = 10e6
T_c = 3800
chamber = H2(P=p_c, T=T_c)


