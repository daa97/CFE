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
        "T_core":3700,
        "T_channel":450,
        "mdot":.108,
        "r1":30e-3,
        "r2":45e-3,
        "r3":49e-3,
        "r5":56e-3,
        "d45":5e-3,
        "d56":3e-3,
        "N":7000,
        "nu_s":0.6961,
        "L_CFE":0.84}

stdlim = [0.5, 2]

# ******************************************
# TODO: if you don't want to plot vs a particular parameter, remove it from `vary`
# TODO: if you want to plot a particular parameter over a range different from others, replace its limits in `vary`
# ******************************************

vary = {"P_core":stdlim,
        "T_channel":stdlim,
        "mdot":stdlim,
        "r5":[0.982,1.5],
        "d56":stdlim,
        "N":stdlim,
        "nu_s":[0.5,1.07],
        "L_CFE":stdlim}

press_vary = {"P_core":stdlim,
            "mdot":stdlim,
            "N":stdlim,
            "L_CFE":stdlim}

pcore = {"P_core":stdlim}
mdot = {"mdot":stdlim}
N = {"N":stdlim}
L = {"L_CFE":stdlim}

labels = {"P_core":"core pressure $P_3$", 
        "T_channel":"channel temperature $T_1$",
        "mdot": "mass flow rate $\dot{m}$",
        "r5":"case radius $r_5$",
        "d56":"outer channel width $(r_6 - r_5)$",
        "N":"CFE rotation rate $\omega$",
        "nu_s":"turbine vel. ratio $\\nu_s$",
        "L_CFE":"CFE length $l$"}

cols = {"P_core":"k",
        "T_channel":"b",
        "r5":"gold",
        "d56":[.09, .75, .19],
        "mdot":"purple",
        "N":"red",
        "nu_s":"gray",
        "L_CFE":"hotpink"}

def run_sweep():
    yvals = dict()
    xvals = dict()

    base_core = H2(P=base["P_core"], T=base["T_core"])      # speed code up by not calculating on every single loop

    for key in vary:                    # iterate through properties we want to vary
        props = base.copy()             # reset all properties to base values
        lim = vary[key]                 # relative property value limits
        n_pts = 50                      # number of x-value points to plot to form a smooth curve
        xvals[key] = np.arange(lim[0], lim[1]+1e-5, np.diff(lim)/n_pts)
        yvals[key] = []
        for x in xvals[key]:
            props[key] = x * base[key]              # adjust single parameter
            L_total = props["L_CFE"] + 0.1          # compute total length
            if key=="P_core" or key=="T_core":      # check if core state needs adjustment
                core = H2(P=props["P_core"], T=props["T_core"])
            else:
                core = base_core
            omega = props["N"] * np.pi/30           # compute omega
            r6 = props["r5"] + props["d56"]         # compute r6


            # ******************************************
            y=1; print("CHANGE THIS LINE!")
            # TODO: ADD YOUR CALCULATIONS FOR THE OUTPUT Y-AXIS VALUE BASED ON INPUTS
            # TODO: SET y=<YOUR Y-VALUE PROPERTY>
            # ******************************************
            
            yvals[key].append(y)


    for key in yvals:       # plot each line
        plt.plot(xvals[key], yvals[key], label=labels[key])

    plt.xlim(0,2.1)
    plt.legend(title="Normalized Parameter", title_fontproperties={"family": "Times New Roman:bold"})
    plt.xlabel("Parameter Normalized to Baseline Configuration")

    # ******************************************
    # TODO: update figure name
    # ******************************************

    plt.savefig("FIGURENAME.svg")
    plt.show()



        