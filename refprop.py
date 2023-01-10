# NOTE: must be run with 32-bit python to work with the free software refpropMINI

import os, numpy as np
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary as RPFL
# ********************************************************************
lib = os.environ['RPPREFIX']
RP = RPFL(lib+"REFPROP.DLL")
RP.SETPATHdll(lib)
SI = RP.GETENUMdll(0,"MASS BASE SI").iEnum
outprops = ["T", "P", "D", "V", "E", "H", "S", "CP", "CP/CV", "W", "M", "VIS", "TCX", "PRANDTL", "G"]

def props(P,T):      # find properties
    r = RP.REFPROPdll("HYDROGEN","TP", ";".join(outprops), SI, 1,0,T,P, [1.])
    assert r.ierr<1, f"Fatal REFPROP error code {RP.ERRMSGdll(r.ierr)} @ P={P}, T={T}"
    o = r.Output[:len(outprops)]           
    return dict(zip(outprops, o))

n=40
trunc = lambda x: round(x,min(-int(np.floor(np.log10(abs(x))))+4,0))
P = np.array([trunc(10**(x/n)) for x in range(2*n,int(7.5*n)+1)]).transpose()
print(max(P)/1e6)
T = np.arange(25, 1505, 5, dtype=int)

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
        if key!="P" and key!="T":
            # loop through every value of the property
            for i in range(len(dictionary["P"])):
                P_idx = np.where(P==np.round(dictionary["P"][i],6))
                T_idx = np.where(T==dictionary["T"][i])
                try:
                    tab[T_idx, P_idx] = dictionary[key][i]
                except:
                    print(key, len(dictionary[key]), len(dictionary["P"]))
                    raise ValueError
            tab = np.concatenate((col1,tab), axis=1)
            tab = np.concatenate((row1,tab))
            tables[key]=tab
    return tables

vals = dict()
for op in outprops:
    vals[op] = []

for Pi in P:
    for Ti in T:
        out = props(P=Pi,T=Ti)
        for k,v in out.items():
            vals[k].append(v)
        if Pi % 100000 == 0 and Ti % 100==0:
            if Ti==100:
                print("*"*50)
            l = [f"{k}={v:.1f}" for k,v in out.items()]
            print(f"|".join(l))



np.savez("refprops.npz", **tabulate(vals, P, T))

