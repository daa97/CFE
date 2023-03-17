import numpy as np
import matplotlib.pyplot as plt
import turbine_design_baines as tdb
import cfe_model as cm
from fluids import FluidsList
import labyrinth_seal_calcs as labby

H2 = FluidsList.H2
air = FluidsList.Air

P1 = np.load("P1.npz")
m_U = np.load("uranium_mass.npz")

opts = {
    "dim" : "Y",
    "prelim" : "y",
    "geom" : "y",
    "stations" : ["y","y"],
    "losses" : "y"
}

nozzle_inputs = {
        "camber angle" : 0.000001,
        "ac" : 0.5,
        "t_2c" : 0.06,#0.1,
        "t_3c" : 0.035,
        "t_maxc" : 0.1,
        "dc" : 0.3,
        "sc" : 0.55,
        "setting angle" : 10/180*np.pi,
        "num_stators" : 20
    }

static_cfe_inputs = {
    "inner_radius" : 0.056, #Channel inner radius [m]
    "outer_radius" : 0.064, #Channel outer radius [m]
    "length" : 0.94, #CFE channel length [m]
    "rpm" : 7000, #CFE inner SiC cylinder revolutions per minute
    "mass_flow" : 0.108, #CFE channel mass flow rate [kg s^-1]
    "uranium_mass":m_U["baseline"],
    "temp" : 450, #[K]
    "press" : 13.763e6, #MPa - Turbine Inlet Pressure
    "fluid" : H2,
    "off_design" : False
} 

test_static_cfe_inputs = {
    "inner_radius" : 0.056, #Channel inner radius [m]
    "outer_radius" : 0.064, #Channel outer radius [m]
    "length" : 0.94, #CFE channel length [m]
    "rpm" : 7000, #CFE inner SiC cylinder revolutions per minute
    "mass_flow" : 0.108, #CFE channel mass flow rate [kg s^-1]
    "uranium_mass":m_U["baseline"],
    "temp" : 450, #[K]
    "press" : 7e6, #MPa - Turbine Inlet Pressure
    "fluid" : H2,
    "inlet_effects" : False
} 

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
    "v_s" : 0.695,
}

test_dynamic_turb_inputs = {
    "PR_ts" : 1.0008,
    "eta_ts" : 0.9,
    "h_0ss" : 0,
    "N_s" : 0,
    "v_s" : 0.708888,
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
    test_turb = tdb.find_turb(init_turb)
    if dict_only:
        test_turb = get_props(test_turb)
        test_cfe = get_props(test_cfe)
    return test_turb

labby_inputs = {
    "fluid" : air,
    "Radial clearance" : 0.000254,# [m]
    "Tooth pitch" : 0.004, # [m]
    "Shaft radius" : 0.112/2,#radius at which the clearance starts
    "Tooth width" : 0.002,
    "Tooth height" : 0.003,
    "Number of teeth" : 7,
    "Inlet pressure" : 15 * 101325,
    "Outlet pressure" : 101325,
    "Inlet temperature" : 287.15,
    "Outlet temperature" : 287.15
}
" 0.0075 < c/s < 0.0375, 0.0075 < w/s < 0.5, 2.67 < w/c < 66.67 and 0.75 < h/s < 4."
if __name__=="__main__":
    # urn_labby = labby.labyrinth_seal(labby_inputs)
    # print("mass flow rate:",urn_labby.mdot,"[kg/s]")
    # print("Pressure profile:",urn_labby.P_profile)
    # print("Density profile:", urn_labby.rho_profile)
    # urn_labby.plot_labby_geom()

    # teeth = [5,7,10,13,15,17,20]
    # PRs = np.linspace(1.01,2.01,30)
    # press = 101325
    # p_given = "outlet"
    # labby.plot_labbys(teeth,press,PRs,air,p_given,labby_inputs)


    # test_cfe = cm.CFE(**static_cfe_inputs)
    # nu_s = np.linspace(0.6,0.71,5)
    # data = tdb.find_lotsa_turbines(nu_s,test_cfe)
    # print(data)

    # Ts = [300,350,400,450,500,550]
    # tdb.plot_lotsa_turbines_temp(Ts)
    # Ps = [10,13.763,17,21,25]
    # tdb.plot_lotsa_turbines_press(Ps)
    # Ms = [0.25,0.5,0.75,0.108,0.125]
    # tdb.plot_lotsa_turbines_mass(Ms)

    test_turb = find_turbine(static_cfe_inputs, dynamic_turb_inputs)
    # # test_turb.make_hub_and_shroud()
    # test_turb.print_turbine(opts)
    # # test_turb.velocity_triangles()

    noz = tdb.nozzle(nozzle_inputs,test_turb)
    noz.export_profile()
    noz.nozzle_eval("y")

    # vals_dict = {
    #     "t_maxc" : [.05,.06,.07,.08,.09,.1,0.125,0.15,0.2],
    #     "t_2c" : [0.025,.03,.04,.04,.05,.06,.075,.1],
    #     "t_3c" : [.012,.025,.03,.035,.04,.045,.05],
    #     "num_stators" : [15,16,17,19,21,23,25,27],
    #     "sc" : [0.5,0.6,0.7,0.75,0.8,0.9]

    # }
    # tdb.stator_sweep_vars(vals_dict,nozzle_inputs,test_turb)