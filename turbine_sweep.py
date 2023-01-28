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

def turb_props(props, P1):
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
    return statics, dynamics

def iter_param(base, key, x, n, parallel=True):
    P1vals = np.load("P1.npz", allow_pickle=True)
    props = base.copy()             # reset all properties to base values
    args = []
    
    for j in range(len(x)):
        props[key] = x[j] * base[key]              # adjust single parameter
        if key in P1vals.keys():
            P1 = P1vals[key][j]
        else:
            P1 = P1vals["baseline"]
        args.append(turb_props(props, P1))

    if parallel:
        pool = Pool()
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
            "r5":stdlim,
            "d56":stdlim,
            "N":stdlim,
            "nu_s":[0.5,1.08],
            "L_CFE":stdlim}

    yvals = dict()
    xvals = dict()
    n_pts = 50
    for key in vary:                    # iterate through properties we want to vary
        lim = vary[key]                 # relative property value limits
        xvals[key] = np.arange(lim[0], lim[1]+1e-5, np.diff(lim)/n_pts)
        yvals[key] = iter_param(base, key, xvals[key], n_pts, parallel=parallel)
    return yvals

if __name__=="__main__":
    turbs = parametric_sweep(parallel=True)
    print("Turbine information collected")

    eta = dict()
    dh = dict()
    radius = dict()
    mass = dict()
    for k,v in turbs.items():
        eta[k]=[]; dh[k]=[]; radius[k]=[]; mass[k]=[];
        for i in v:
            eta[k].append(i["eta_ts_loss"])
            dh[k].append(i["deltah_0"])
            radius[k].append(i["r_4"])
    
    np.savez("eta_sweep_vs.npz", **eta)
    np.savez("dh_sweep_vs.npz", **dh)
    np.savez("radius_sweep_vs.npz", **radius)

    print("Done!")

    


            

