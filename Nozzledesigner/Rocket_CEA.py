# **********Step 1: Create CEA input.inp files***********
import os
from random import randint as rint

def create_CEA_infile(p_c,T_c, area_rat,supersonic,debug):
    os.chdir(r"C:\Users\tdham\OneDrive - Georgia Institute of Technology\Years\Third\Fall\NASA\CFE-Thermal-main" )
    dir = "RCEAexec\\CEAexec-win\\"
    if supersonic == False:
        in_template ="rcea_sub.txt"
    else:
        in_template ="rcea_sup.txt"
    with open(dir+in_template, mode='r') as f_in_template:
        contents = f_in_template.read()
        if supersonic == False:
            contents = contents.replace("<SUBSONIC_AREA_RATIO>", str(area_rat))
        else:
            contents = contents.replace("<SUPERSONIC_AREA_RATIO>", str(area_rat))
        contents = contents.replace("<CHAMBER_STAGNATION_TEMPERATURE>", str(T_c))
        contents = contents.replace("<CHAMBER_STAGNATION_PRESSURE>", str(p_c))
        writename = f"Temp_RCEA_in{rint(0,1000)}"
    with open(dir+writename+'.inp', mode='w') as f_in:
        f_in.write(contents)

        if debug:
            with open(dir + writename + '.txt', mode='w') as f_in_txt:
                f_in_txt.write(contents)
            with open(dir + writename + '.txt', mode='r') as f_in_txt:
                print(f_in_txt.read())

    return writename

# Step 2: Run CEA and create CEA output file

# Step 3: Parse Output File; store in arrays

