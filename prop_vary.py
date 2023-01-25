from fluids import H2
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
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
        "nu_s":0.691,
        "L_CFE":.84,
        "T_core":3700}



stdlim = [0.5, 2]


# TODO: if you don't want to plot vs a particular parameter, remove it from `vary`
# TODO: if you want to plot a particular parameter over a range different from others, replace its limits in `vary`
vary = {"P_core":stdlim,
        "T_channel":stdlim,
        "r5":stdlim,
        "d56":[0.125, 2],
        "N":stdlim,
        "nu_s":stdlim,
        "L_CFE":stdlim}

labels = {"P_core":"core pressure $P_3$", 
        "T_channel":"channel temperature $T_1$",
        "r5":"case radius $r_5$",
        "d56":"outer channel width $(r_6 - r_5)$",
        "N":"CFE rotation rate $\omega$",
        "nu_s":"turbine sp. speed $\\nu_s$",
        "L_CFE":"CFE length $l$"}

yvals = dict()
base_core = H2(P=base["P_core"], T=base["T_core"])
for key in vary:
    props = base.copy()             # reset all properties to base values
    lim = vary[key]                 # relative property value limits
    n_pts = 50                      
    points = np.arange(lim[0], lim[1]+1e-5, np.diff(lim)/n_pts)
    yvals[key] = []
    for pt in points:
        props[key] = pt * base[key]
        L_total = props["L_CFE"] + 0.1
        if key=="P_core" or key=="T_core":
            core = H2(P=props["P_core"], T=props["T_core"])
        else:
            core = base_core
        omega = props["N"] * np.pi/30
        r6 = props["r5"] + props["d56"]


        # ******************************************
        y=1
        # TODO: ADD YOUR CALCULATIONS FOR THE OUTPUT Y-AXIS VALUE BASED ON INPUTS
        # TODO: SET y=<YOUR Y-VALUE PROPERTY>

        # ******************************************
        
        yvals[key].append(y)


for key in yvals:
    plt.plot(yvals[key], label=labels[key])

plt.legend()    
plt.xlim(0,2.2)
plt.show()



        