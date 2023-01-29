import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import testforhs as tfhs
from multiprocessing import Pool
from os import cpu_count


mpl.rc('font', family='Times New Roman',size="10")
mpl.rc('figure', figsize=(4.8,3.6))
mpl.rc('savefig', dpi=800)
mpl.rc('lines', linewidth=1.2)
mpl.rc('axes', grid=True)
mpl.rc('grid', linewidth=0.25)
mpl.rc('mathtext', fontset="dejavuserif")
mpl.rc('xtick.minor', visible=True, size=1.5, width=0.5)
mpl.rc('ytick.minor', visible=True, size=1.5, width=0.5)
plt.rcParams['figure.constrained_layout.use'] =  True

base = {"P_core":10e6,
        "T_channel":450,
        "r5":56e-3,
        "d56":8e-3,
        "N":7000,
        "nu_s":0.693,
        "L_CFE":.84,
        "T_core":3700}

stdlim = [0.5, 2]

vary = {"P_core":stdlim,
        "T_channel":stdlim,
        "r5":[0.93,1.5],
        "d56":stdlim,
        "N":stdlim,
        "nu_s":[0.5,1.08],
        "L_CFE":stdlim}


P1vals = np.load("P1.npz", allow_pickle=True)
mvals = np.load("uranium_mass.npz", allow_pickle=True)


def turb_props(props, P1, m):
    """Takes a property dictionary built for performing a generic parametric
    sweep and turns it into a set of inputs for the find_turb function."""
    L_total = props["L_CFE"] + 0.1          # compute total length
    r6 = props["r5"] + props["d56"]         # compute r6
    statics = tfhs.static_cfe_inputs.copy()
    dynamics = tfhs.dynamic_turb_inputs.copy()
    statics["inner_radius"] = props["r5"]
    statics["outer_radius"] = r6
    statics["rpm"] = props["N"]
    statics["temp"] = props["T_channel"]
    statics["length"] = L_total
    dynamics["v_s"] = props["nu_s"]
    statics["press"] = P1 / 1e6
    return statics, dynamics, True

def iter_param(base, key, x, n, parallel=True):
    args = []
    props = base.copy()             # reset all properties to base values

    for j in range(len(x)):
        props[key] = x[j] * base[key]              # adjust single parameter
        if key in P1vals.keys():
            P1 = P1vals[key][j]
            m = mvals[key][j]
        else:
            P1 = P1vals["baseline"]
            m = mvals["baseline"]
        args.append(turb_props(props, P1, m))

    if parallel:
        pool = Pool(14)
        out = pool.starmap(tfhs.find_turbine, args)
        pool.close()
        pool.join()
    else:
        out = []
        for j in range(len(args)):
            print("*"*10, f"{x[j]:.2f} [{j}/{n}]", "*"*10)
            out.append(tfhs.find_turbine(*args[j]))

    return out

def parametric_sweep(parallel=True):
    yvals = dict()
    xvals = dict()
    n_pts = 50
    for key in vary:                    # iterate through properties we want to vary
        lim = vary[key]                 # relative property value limits
        xvals[key] = np.arange(lim[0], lim[1]+1e-5, np.diff(lim)/n_pts)
        yvals[key] = iter_param(base, key, xvals[key], n_pts, parallel=parallel)
    return yvals

if __name__=="__main__":
    turb_cfes = parametric_sweep(parallel=True)
    print("Turbine information collected")
    turb_cfes["baseline"] = [tfhs.find_turbine(*turb_props(base, P1vals["baseline"], mvals["baseline"]))]
    
    savenames = ['eta', 'W_bear', 'W_visc', 'W', 'radius', 'mass']
    outs = dict()
    for s in savenames:
        outs[s] = dict()
    for k,v in turb_cfes.items():
        for s in outs:
            outs[s][k] = []
        for turb, cfe in v:
            outs['eta'][k].append(turb["eta_ts_loss"])
            outs['radius'][k].append(turb["r_4"])
            outs['W'][k].append(cfe["work_rate"])
            outs['W_bear'][k].append(cfe["M_bearing"] * cfe["omega"])
            outs['W_visc'][k].append(cfe["M_visc"] * cfe["omega"])
            outs['mass'][k].append(cfe["mass"])

    for s in outs:
        np.savez(f"turbine_sweep/{s}.npz", **outs[s])
    
    print("Done!")

    


            

