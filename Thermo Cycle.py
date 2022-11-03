from cProfile import label
from fluids import *

import matplotlib.pyplot as plt
def moment(mass, R1, R2):
    return (1/2) * mass * (R1**2 + R2**2)

def mass(R1, R2, length, density):
    return np.abs(np.pi * (R2**2 - R1**2) * length) * density


# ***** Dimensions and Power Calc Parameters ************
l = 0.84
r1 = 0.03
r2 = 0.045
rho_u = 17100
void_frac = 0.5
omega = 7000 * (2 * np.pi / 60)
t_init = 60*30 # seconds for spin-up
cfe_flow = 0.1
n = 19 # Number of CFEs

# ********** Power/Force Calculations *******************
rho_uh = (1 - void_frac) * rho_u
mass_u = mass(r1, r2, l, rho_uh)
moment_u = moment(mass_u, r1, r2)
ke_cfe = (1/2) * moment_u * omega**2
power_per_bearing = 50
bearing_loss = power_per_bearing * 3
cfe_power =(ke_cfe / t_init) + bearing_loss
F_cent_ur = 2/3 * omega**2 * (r2**3 - r1**3) * np.pi * rho_uh * l
A_surf_CFE = 2 * np.pi * r2 * l

P_cent = F_cent_ur/(A_surf_CFE)

# *********** PM Pressure Drop ********************
P_chamber = 10000                                                                   # Chamber pressure (kPa) (why is this different from p_core?)
T_pm_2 = 1000                                                                       # Temperature at PM (again different from T_pm)
R_H2 = 4.124                                                                        # Hydrogen gas constant
rho_H2 = P_chamber / (R_H2 * T_pm_2)                                                # Density of H2
mu_H2_PM = -(0.00000000000144)*T_pm_2 + (0.0000000169)*T_pm_2 + 0.00000464          # Viscosity at PM (Pa-s)
r_outer_PM = 0.049                                                                  # Outer radius of PM
r_inner_PM = 0.0455                                                                 # Inner radius of PM
L_PM = r_outer_PM - r_inner_PM                                                      # Thickness of PM
D_load_PM = 2 * P_chamber * L_PM                                                    # Distributed load across length of PM (N/m)
m_flow_system = 2.3                                                                 # Mass flow of system
cfe_flow2 = m_flow_system / n                                                       # Mass flow per CFE (not equal to cfe_flow above?)
porosity_PM = 0.3                                                                   # PM volume fraction
cfe_length = 0.81                                                                   # CFE length
q_H2 = cfe_flow2 / (rho_H2 * porosity_PM * cfe_length * 2 * np.pi * r_outer_PM)     # Darcian Velocity (specific discharge)
k1 = 3.05322463490072E-09                                                                                   # Permeability, k1
k2 = 0.000421605940239228                                                                                   # Permeability, k2
P_PM_inlet = np.sqrt(D_load_PM * ((mu_H2_PM/1000) * (q_H2/k1) + rho_H2 * (q_H2)**2 / k2) + P_chamber**2)    # PM inlet pressure
dp_pm2 = P_chamber - P_PM_inlet                                                                             # Pressure drop across PM

# *********** Thermo Fluid Parameters **************
T_tank = 30             # starting temperature
P_tank = 812930.475     # starting pressure (vapor pressure of H2 @ 30 K)
T_core = 3700           # final temperature inside CFE
p_core = 5e6            # final pressure inside CFE
T_pm = 500              # temperature at porous media
dp_uran = -P_cent       # pressure change across uranium due to centripetal pressure
# dp_pm = -9e3            # pressure change across porous media
dp_pm = dp_pm2            # pressure change across porous media
eta_turbopump = 1       # turbopump efficiency
eta_cfeturb = 1         # cfe turbine efficiency

dh_cfeturb = - cfe_power / cfe_flow 
dh_turbo = -331706.6         # starting value for turbopump work, will be iteratively solved later
y_turb = 0.5              # fraction going into turbopump
q_regen1 = 0            # turbopump bypass heating states 3->6
q_regen2 = 0            # CFE entrance heating states 5->8
y_throt = 1-y_turb      # fraction which bypasses via throttle

H2 = Fluid("Hydrogen", prop_files)
start = H2(P=P_tank, T=T_tank)                          # state 1

core = H2(T=T_core, P=p_core)                           # state 11
pm = H2(T=T_pm, P=core.P - dp_uran)                     # state 10
cfe = H2(h=pm.h, P=pm.P - dp_pm)                        # state 9
hot_isent = H2(h=cfe.h-dh_cfeturb, s=cfe.s)             # state 8s
hot = H2(h=cfe.h-dh_cfeturb*eta_cfeturb, P=hot_isent.P) # state 8
mix = H2(h=hot.h-q_regen2, P=hot.P)                     # state 5


# TODO: **IMPORTANT** Account for TP bypass split and remixing
# TODO: Update starting temp/pressure to vapor pressure
# TODO: Update PM dP and uranium dP to reflect accurate values
# TODO: Account for inefficient TP (both halves separately?) and check CFE efficiency calcs
# TODO: Maybe add parametric sweeps for TP bypass %, as well as 2nd and 3rd regen stages


while True:
    turbo = H2(h=mix.h + y_throt*dh_turbo, P=mix.P)     # state 4
    throt = H2(h=mix.h - y_turb*dh_turbo, P=mix.P)      # state 7
    regen = H2(h=turbo.h-dh_turbo, s=turbo.s)           # state 3
    bypass = H2(h=throt.h, P=regen.P)                   # state 6
    pumped = H2(P=regen.P, s=start.s)                   # state 2

    dh_pump_goal = - dh_turbo / y_turb / eta_turbopump  # expected pump work based upon turbine
    dh_pump = pumped.h - start.h                        # empirical pump work
    q_nozzle = regen.h - pumped.h                       # required regenerative cooling
    err = dh_pump - dh_pump_goal                        # difference of expcted and actual pump work found
    print("-"*40)
    print("Converging ---- enthalpy error:", err)
    print("Turbopump enthalpy change:", dh_pump)
    print("Regen nozzle cooling:", q_nozzle)
    print("Max Pressure:", pumped.P)
    if abs(err) < 1:
        print("Converged!")
        break
    dh_turbo -= err * y_turb

bypass = H2(h=regen.h+q_regen1, P=regen.P)              # state 6
throt = H2(h=bypass.h, P=mix.P)                        # state 5

space = H2(s=core.s, T=300)
print(space.P)

#states = [start, pumped, regen, turbo, mix, bypass, throt, hot, cfe, pm, core]
states = [start, pumped, regen, turbo, throt, mix, cfe, pm, core] # modified to not plot duplicates
flow1 = [start, pumped, regen, bypass, throt, mix, hot, cfe, pm, core, space]
flow2 = [regen, turbo, mix]
n = 2
f1 = []
for i in range(len(flow1)-1):
    p1 = flow1[i]
    f1.append(p1)
    p2 = flow1[i+1]
    for j in range(1,n):
        if p1.s == p2.s:
            pass
            #f1.append(H2(s=p1.s, h=(p2.h-p1.h)*j/n + p1.h))
        elif p1.h == p2.h:
            pass
            #f1.append(H2(h=p1.h, s=(p2.s-p1.s)*j/n + p1.s))
        else:
            f1.append(H2(h=(p2.h-p1.h)*j/n + p1.h, P=(p2.P-p1.P)*j/n + p1.P))
f1.append(flow1[-1])


fig, axes = plt.subplots(2,1)
for ax in axes:
    ax.plot([point.s/1e3 for point in f1], [point.T for point in f1], 'k')
    ax.plot([point.s/1e3 for point in flow2], [point.T for point in flow2], 'k')
    ax.grid(True)
    styles = '.^v+x'
    for i, point in enumerate(states):
        ax.plot(point.s/1e3, point.T, styles[i%5], label=i+1)
    ax.legend()
    #plt.ylabel("Temperature (K)")
    #plt.xlabel("Entropy (kJ/kg K)")

fig.show()

for point in states:
    print(point.T, "|", point.s)




