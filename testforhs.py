import numpy as np
import matplotlib.pyplot as plt
import turbine_design_baines as tdb
import cfe_model as cm
import prop_vary as pv
from fluids import FluidsList
# etas = [0.78, 0.80, 0.82, 0.84, 0.85 0.86]
# nus = range(0.15,1,0.01)
# for eta in etas:
#     for nu in nus:
#         psi = eta/(2 * nu**2)

P1 = np.load("P1.npz")
m_U = np.load("uranium_mass.npz")

opts = {
    "dim" : "Y",
    "prelim" : "y",
    "geom" : "y",
    "stations" : ["Y","y"],
    "losses" : "Y"
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

# static_cfe_inputs = {
#     "inner_radius" : 0.056, #Channel inner radius [m]
#     "outer_radius" : 0.064, #Channel outer radius [m]
#     "length" : 0.94, #CFE channel length [m]
#     "rpm" : 7000, #CFE inner SiC cylinder revolutions per minute
#     "mass_flow" : 0.108, #CFE channel mass flow rate [kg s^-1]
#     "uranium_mass":m_U["baseline"],
#     "temp" : 450, #[K]
#     "press" : P1["baseline"], #Pa - Turbine Inlet Pressure
#     "fluid" : FluidsList.Air
# }

off_design_CFE = {
    "inner_radius" : 0.056, #Channel inner radius [m]
    "outer_radius" : 0.0635, #Channel outer radius [m]
    "length" : 0.50, #CFE channel length [m]
    "rpm" : 1000, #CFE inner SiC cylinder revolutions per minute
    "mass_flow" : 0.003, #CFE channel mass flow rate [kg s^-1]
    "uranium_mass":0,
    "temp" : 300, #[K]
    "press" : 101.325e3, #Pa - Turbine Inlet Pressure
    "fluid" : FluidsList.Air,
    "inlet_effects" : True
}

dynamic_turb_inputs = {
    "PR_ts" : 1.0008,
    "eta_ts" : 0.9,
    "h_0ss" : 0,
    "N_s" : 0,
    "v_s" : 0.695
}

def get_props(turb):
    copy = dict()
    for k,v in turb.__dict__.items():
        if isinstance(v, (str, int, float, np.ndarray)):
            copy[k] = v
    return copy

def find_turbine(static_inputs, dynamic_inputs, dict_only=False):
    test_cfe = cm.CFE(**static_inputs)
    init_turb = tdb.turbine(test_cfe.static_turb_inputs,dynamic_inputs,1)
    test_turb = tdb.find_turb(test_cfe,init_turb)
    if dict_only:
        test_turb = get_props(test_turb)
        test_cfe = get_props(test_cfe)
    return test_turb, test_cfe

if __name__=="__main__":
    
    test_turb, test_cfe = find_turbine(static_inputs=off_design_CFE, dynamic_inputs=dynamic_turb_inputs, dict_only=True)
    print(test_cfe)