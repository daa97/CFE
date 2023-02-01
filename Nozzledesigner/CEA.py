import os
import subprocess as sp
from config.definitions import ROOT_DIR

def runCEA(file):
    os.chdir(os.path.join(ROOT_DIR, 'RCEAexec', 'CEAexec-win'))
    sp_status = sp.run(["echo", f"{file}", "|", "FCEA2.exe"], stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
    return sp_status