import numpy as np
import matplotlib.pyplot as plt
import turbine_design_baines as tdb
# etas = [0.78, 0.80, 0.82, 0.84, 0.85 0.86]
# nus = range(0.15,1,0.01)
# for eta in etas:
#     for nu in nus:
#         psi = eta/(2 * nu**2)
static_cfe_inputs = {
    "inner_radius" : 0.056, #Channel inner radius [m]
    "outer_radius" : 0.064, #Channel outer radius [m]
    "length" : 0.94, #CFE channel length [m]
    "rpm" : 7000, #CFE inner SiC cylinder revolutions per minute
    "mass_flow" : 0.108, #CFE channel mass flow rate [kg s^-1]

    "temp" : 450, #[K]
    "press" : 13.1135 #MPa - Turbine Inlet Pressure
} 

dynamic_turb_inputs = {
    "PR_ts" : 1.0008,
    "eta_ts" : 0.9,
    "h_0ss" : 0,
    "N_s" : 0,
    "v_s" : 0.693
}
opts = {
    "dim" : "Y",
    "prelim" : "y",
    "geom" : "y",
    "stations" : ["Y","y"]
}

nozzle_inputs = {
        "radius ratio" : 1.1,
        "camber angle" : np.pi/6,
        "ac" : 0.25,
        "t_2c" : 0.025,
        "t_3c" : 0.012,
        "t_maxc" : 0.06,
        "dc" : 0.4,
        "sc" : 0.75,
        "setting angle" : 10/180*np.pi
    }

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

static_cfe_inputs = {
    "inner_radius" : 0.056, #Channel inner radius [m]
    "outer_radius" : 0.064, #Channel outer radius [m]
    "length" : 0.94, #CFE channel length [m]
    "rpm" : 7000, #CFE inner SiC cylinder revolutions per minute
    "mass_flow" : 0.108, #CFE channel mass flow rate [kg s^-1]

    "temp" : 450, #[K]
    "press" : 13.1135 #MPa - Turbine Inlet Pressure
} 

dynamic_turb_inputs = {
    "PR_ts" : 1.0008,
    "eta_ts" : 0.9,
    "h_0ss" : 0,
    "N_s" : 0,
    "v_s" : 0.693
}
def find_turbine(static_inputs, dynamic_inputs):
    test_cfe = tdb.CFE(static_inputs,dynamic_inputs,1)
    init_turb = tdb.turbine(test_cfe,test_cfe.static_turb_inputs,dynamic_inputs,1)
    test_turb = tdb.find_turb(test_cfe,init_turb)
    return test_turb
def parametric_sweep():
    base = {"P_core":10e6,
            "T_channel":450,
            "r5":56e-3,
            "d56":8e-3,
            "N":7000,
            "nu_s":0.691,
            "L_CFE":.84,
            "T_core":3700}

    stdlim = [0.5, 2]

    # ******************************************
    # TODO: if you don't want to plot vs a particular parameter, remove it from `vary`
    # TODO: if you want to plot a particular parameter over a range different from others, replace its limits in `vary`
    # ******************************************

    vary = {"P_core":stdlim,
            "T_channel":stdlim,
            "r5":stdlim,
            "d56":stdlim,
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
    xvals = dict()

    base_core = H2(P=base["P_core"], T=base["T_core"])      # speed code up by not calculating on every single loop
    for key in vary:                    # iterate through properties we want to vary
        props = base.copy()             # reset all properties to base values
        lim = vary[key]                 # relative property value limits
        n_pts = 50                      # number of x-value points to plot to form a smooth curve
        xvals[key] = np.arange(lim[0], lim[1]+1e-5, np.diff(lim)/n_pts)
        
        print(xvals)
        yvals[key] = []
        P1vals = np.load("P1.npz", allow_pickle=True)
        for j in range(len(xvals[key])):
            props[key] = xvals[key][j] * base[key]              # adjust single parameter
            L_total = props["L_CFE"] + 0.1          # compute total length
            if key=="P_core" or key=="T_core":      # check if core state needs adjustment
                core = H2(P=props["P_core"], T=props["T_core"])
            else:
                core = base_core
            omega = props["N"] * np.pi/30           # compute omega
            r6 = props["r5"] + props["d56"]         # compute r6
            if key in P1vals.keys():
                P1 = P1vals[key][j]
            else:
                P1 = P1vals["baseline"]
            statics = static_cfe_inputs.copy()
            dynamics = dynamic_turb_inputs.copy()
            statics["inner_radius"] = props["r5"]
            statics["outer_radius"] = r6
            statics["rpm"] = props["N"]
            statics["temp"] = props["T_channel"]
            statics["length"] = L_total
            dynamics["v_s"] = props["nu_s"]
            statics["press"] = P1
            turb = find_turbine(statics, dynamics)
            yvals[key].append(turb)
    return yvals


turbs = parametric_sweep()


test_turb = find_turbine(static_inputs=static_cfe_inputs, dynamic_inputs=dynamic_turb_inputs)
test_turb.make_hub_and_shroud()
# test_turb.print_turbine(opts)

# noz = tdb.nozzle(nozzle_inputs,test_turb)

# noz.create_cascade()
# noz.find_setting_angle()
# noz_out = noz.calc_naca_profile()
# nozx = noz_out[0]
# nozy = noz_out[1]
# nozchi = noz_out[2]
# nozsuc = noz_out[3]
# nozpres = noz_out[4]
# r_2 = noz.t_2c/2 * noz.c
# r_3 = noz.t_3c/2 * noz.c
# thetas = np.linspace(0, 2 * np.pi, num = 100)
# le = np.zeros((100,2))
# te = np.zeros((100,2))
# for i,theta in enumerate(thetas):
#     le[i,0] = r_2 * (np.cos(theta))
#     le[i,1] = r_2 * (np.sin(theta))
#     te[i,0] = r_3 * (np.cos(theta)) + nozx[-1]
#     te[i,1] = r_3 * (np.sin(theta)) + nozy[-1]
# plt.plot(le[:,0],le[:,1])
# plt.plot(te[:,0],te[:,1])
# plt.plot(nozx,nozy,marker = ".")
# plt.plot(nozpres[:,0],nozpres[:,1],marker = ".")
# plt.plot(nozsuc[:,0],nozsuc[:,1],marker = ".")
# # plt.xlim([-0.01,1.01])
# # plt.ylim([-0.1,0.1])
# plt.show()
# test_cfe = tdb.CFE(static_cfe_inputs,dynamic_turb_inputs,1)

# init_turb = tdb.turbine(test_cfe,test_cfe.static_turb_inputs,dynamic_turb_inputs,1)

# # test_turb = tdb.find_turb(test_cfe,init_turb)
# # test_turb.print_turbine(opts)
# init_turb.make_hub_and_shroud()