# **********Step 2: Run CEA to produce output.out files***********
import os
import subprocess as sp

def runCEA(file):

    if not ("CEAexec-win" in os.getcwd()):
        os.chdir(r"C:\Users/\tdham\OneDrive - Georgia Institute of Technology\Years\Third\Fall\NASA\CFE-Thermal-main\CEAexec\CEAexec-win")

    sp_status = sp.run(["echo", f"{file}", "|", "FCEA2.exe"], stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
    return sp_status