import numpy as np
import os
import subprocess as sp
import csv

lib = os.environ['CEA']

def match(props, conv, line):
    for key, value in props.items():
        if key in line and len(line)>16:
            numtext = line[16:]
            columns=[]
            for i in range(0,len(numtext)-8,9):
                columns.append(numtext[i:i+9])
            for property_text in columns:
                if property_text[-2]=='-' or property_text[-2]==' ':
                    try:
                        propval = float(property_text[:-2])*10**float(property_text[-2:])
                    except ValueError:
                        if "****" in property_text:
                            propval = 1.0*10**float(property_text[-2:])
                        else:
                            raise ValueError
                else:
                    propval = float(property_text)
                value.append(propval*conv[key])
    return props


def write_inp(P: list, T: list, batch:int=15, template:str="template.inp", prefix:str="PY_") -> list:
    fnames = []
    with open(lib+"\\"+template, mode='r') as f:
        blank = f.read()
    for i in range(int(np.ceil(len(P)/batch))):
        for j in range(int(np.ceil(len(T)/batch))):
            Pi = P[i*batch:min(len(P),(i+1)*batch)]
            Tj = T[j*batch:min(len(T),(j+1)*batch)]
            Pstr = ','.join([str(pi) for pi in np.round(Pi/1e5,6)])
            Tstr = ','.join([str(ti) for ti in Tj.tolist()])
            contents = blank.replace("<PYPRES>", Pstr).replace("<PYTEMP>", Tstr)
            writename = f"{prefix}P{i:02d}T{j:02d}"
            fnames.append(writename)
            with open(lib+"\\"+writename+'.inp', mode='w') as fwrite:
                fwrite.write(contents)
    return fnames

def run_cea(fnames):
    os.chdir(lib)
    for file in fnames:
        sub = sp.run(["echo", f"{file}", "|", "FCEA2.exe"], stdout=sp.PIPE, stderr=sp.PIPE, shell=True)

def read_out(fnames):
    thermo = {"P, BAR":[], "T, K":[], "RHO, KG/CU M":[], "H, KJ/KG":[], "U, KJ/KG":[], 
            "G, KJ/KG":[], "S, KJ/(KG)(K)":[], "M, (1/n)":[], "GAMMAs":[], "SON VEL,M/SEC":[], "VISC,MILLIPOISE":[]}

    thermo_conv = {"P, BAR":1e5, "T, K":1, "RHO, KG/CU M":1, "H, KJ/KG":1e3, "U, KJ/KG":1e3, 
            "G, KJ/KG":1e3, "S, KJ/(KG)(K)":1e3, "M, (1/n)":1e-3, "GAMMAs":1, "SON VEL,M/SEC":1, "VISC,MILLIPOISE":1e-4}
    frozenheader = "WITH FROZEN REACTIONS"
    frozenend = "THERMODYNAMIC PROPERTIES"
    frozen = False
    transp = {"Cp, KJ/(KG)(K)":[], "CONDUCTIVITY   ":[],"PRANDTL NUMBER":[]}
    transp_conv = {"Cp, KJ/(KG)(K)":1e3, "CONDUCTIVITY   ":1e-1,"PRANDTL NUMBER":1}

    # iterate over output files
    io=0
    for f in fnames:
        # open each file
        with open(f+".OUT") as file:
            lines = file.readlines()
        # iterate over each line of the file
        for num in range(len(lines)):
            textline = lines[num]
            if frozenheader in textline:
                frozen = True
            elif frozenend in textline:
                frozen = False
            # go through properties
            thermo = match(thermo, thermo_conv, textline)
            if frozen:
                transp = match(transp, transp_conv, textline)

    thermo.update(transp)
    return thermo

def tabulate(dictionary, P, T):
    # set up first row and first column containing pressure and temperature values
    col1 = np.zeros((len(T),1))
    col1[:,0]=T
    row1 = np.zeros((1,len(P)+1))
    row1[0,1:]=P
    tables = dict()
    # loop through each property
    for key in dictionary.keys():
        # create an empty array for property
        tab = np.zeros((len(T), len(P)))
        tab[:] = np.nan
        # do not generate pressure and temperature tables as these serve as our axes
        if key!="P, BAR" and key!="T, K":
            # loop through every value of the property
            for i in range(len(dictionary["P, BAR"])):
                P_idx = np.where(P==np.round(dictionary["P, BAR"][i],2))
                T_idx = np.where(T==dictionary["T, K"][i])
                try:
                    tab[T_idx, P_idx] = dictionary[key][i]
                except ValueError:
                    print(key, len(dictionary[key]), len(dictionary["P, BAR"]))
                    raise ValueError
            tab = np.concatenate((col1,tab), axis=1)
            tab = np.concatenate((row1,tab))
            
            name = key.split(",")[0].strip().replace(" ", "_")
            filt = filter(lambda s: str.isalpha(s) or s=="_", name)
            newkey = ''.join(filt)
            tables[newkey]=tab
    return tables

def write_csvs(tables):
    for key, tab in tables.items():
        csv_fname = key + ".csv"
        with open(csv_fname, mode='w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerows(tab)