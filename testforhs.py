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

def get_props(turb):
    copy = dict()
    for k,v in turb.__dict__.items():
        if isinstance(v, (str, int, float, np.ndarray)):
            copy[k] = v
    return copy

def find_turbine(static_inputs, dynamic_inputs):
    test_cfe = tdb.CFE(static_inputs,dynamic_inputs,1)
    init_turb = tdb.turbine(test_cfe,test_cfe.static_turb_inputs,dynamic_inputs,1)
    test_turb = tdb.find_turb(test_cfe,init_turb)
    return get_props(test_turb)

if __name__=="__main__":

    test_turb = find_turbine(static_inputs=static_cfe_inputs, dynamic_inputs=dynamic_turb_inputs)
    copy = dict()
    for k,v in test_turb.__dict__.items():
        if isinstance(v, (str, int, float, np.ndarray)):
            copy[k] = v
    
    for k,v in copy.items():
        print(f"'{k}' : {v}")
    # test_turb.make_hub_and_shroud()

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