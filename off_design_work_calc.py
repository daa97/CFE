import numpy as np
import matplotlib.pyplot as plt
import turbine_design_baines as tdb
import cfe_model as cm
import prop_vary as pv
from fluids import FluidsList
import csv
import off_design_pressure as odp

Air = FluidsList.Air
def turb_outlet_pressure(core_pressure, mdot, R2, R3, L, T_PM):
    "Returns turbine outlet pressure P1 and SiC inner wall pressure P2"
    P2 = core_pressure
    #print("P2", P2)
    PMstate = Air(T=T_PM, P=P2) # Fluid state
    Q= mdot / PMstate.rho # Volumetric flow rate [m^3/s]
    A = 2*np.pi*R3*L # Cross sectional area [m^2]
    vs = Q / A # Darcian velocity [m/s]
    k1 = 5.69771546674471E-12 # for 18% open porosity, 35% total
    k2 = np.exp(-1.71588/(k1**0.08093)) # Correlation from Dey et al.
    P1 = np.sqrt((PMstate.mu*vs/k1 + (PMstate.rho*vs**2)/k2) * 2*P2*(R3-R2) + P2**2) # Inlet pressure to PM [Pa]
    PM_dP =  P1 - P2 # Pressure drop across PM [Pa]
    return P1, P2

P1 = np.load("P1.npz")
m_U = np.load("uranium_mass.npz")

opts = {
    "dim" : "Y",
    "prelim" : "y",
    "geom" : "y",
    "stations" : ["Y","y"],
    "losses" : "Y"
}
mdot = .014
off_design_CFE = {
    "inner_radius" : 0.056, #Channel inner radius [m]
    "outer_radius" : 0.0635, #Channel outer radius [m]
    "length" : 0.50, #CFE channel length [m]
    "rpm" : 1000, #CFE inner SiC cylinder revolutions per minute
    "mass_flow" : mdot, #CFE channel mass flow rate [kg s^-1]
    "uranium_mass":0,
    "temp" : 300, #[K]
    "press" : odp.outlet_pressure(mdot), #Pa - Turbine Inlet Pressure
    "fluid" : FluidsList.Air,
    "off_design" : True
}


def get_props(turb):
    copy = dict()
    for k,v in turb.__dict__.items():
        if isinstance(v, (str, int, float, np.ndarray)):
            copy[k] = v
    return copy

if __name__=="__main__":
    x = []; y = []; y1 = []; y2 = []; y3 = []
    for i in range(200):
        x.append(i*1e-4)
        off_design_CFE["mass_flow"] = x[-1]
        off_design_CFE["press"] = odp.outlet_pressure(x[-1])
        C = cm.CFE(**off_design_CFE)
        y.append(C.work_rate)
        y1.append(C.M_bearing*C.omega)
        y2.append(C.M_visc*C.omega)
        y3.append(C.M_inlet*C.omega)

    plt.plot(x, y, label="total")
    plt.plot(x, y1, "--", label="Bearing")
    plt.plot(x, y2, "--", label="TTCF")
    plt.plot(x, y3, "--", label="Inlet")
    with open("off_design_work.csv", mode='w', newline="\n") as f:
        w = csv.writer(f)
        w.writerow(["mdot (kg/s)", "power (W)"])
        w.writerows(zip(x,y))
    plt.legend()
    plt.title("Work vs. Mass flow")
    plt.show()