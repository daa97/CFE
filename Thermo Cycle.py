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

# ********** Power Calculations *******************
rho_uh = (1 - void_frac) * rho_u
mass_u = mass(r1, r2, l, rho_uh)
moment_u = moment(mass_u, r1, r2)
ke_cfe = (1/2) * moment_u * omega**2
power_per_bearing = 50
bearing_loss = power_per_bearing * 3
cfe_power =(ke_cfe / t_init) + bearing_loss


# *********** Thermo Fluid Parameters **************
P_tank = 1.5e5          # starting pressure
T_tank = 40             # starting temperature
p_core = 5e6            # final pressure
T_core = 3700           # final temperature
T_pm = 500              # temperature at porous media
dp_uran = -2.2e6        # pressure change across uranium due to centripetal pressure
dp_pm = -9e3            # pressure change across porous media
eta_turbopump = 1       # turbopump efficiency
eta_cfeturb = 1         # cfe turbine efficiency

dh_cfeturb = - cfe_power / cfe_flow 
dh_turbo = -6e6         # starting value for turbopump work, will be iteratively solved later
y_turb = 1              # fraction going into turbopump
q_regen1 = 0            # turbopump bypass heating states 3->6
q_regen2 = 0            # CFE entrance heating states 5->8
y_throt = 1-y_turb      # fraction which bypasses via throttle

H2 = Fluid("Hydrogen", prop_files)

core = H2(T=T_core, P=p_core)                           # state 11
pm = H2(T=T_pm, P=core.P - dp_uran)                     # state 10
cfe = H2(h=pm.h, P=pm.P - dp_pm)                        # state 9
hot_isent = H2(h=cfe.h-dh_cfeturb, s=cfe.s)             # state 8s
hot = H2(h=cfe.h-dh_cfeturb*eta_cfeturb, P=hot_isent.P) # state 8
mix = H2(h=hot.h-q_regen2, P=hot.P)                     # state 8
turbo = H2(h=mix.h, P=mix.P)                            # state 4
start = H2(P=P_tank, T=T_tank)                          


# TODO: **IMPORTANT** Account for TP bypass split and remixing
# TODO: Update starting temp/pressure to vapor pressure
# TODO: Update PM dP and uranium dP to reflect accurate values
# TODO: Account for inefficient TP (both halves separately?) and check CFE efficiency calcs
# TODO: Maybe add parametric sweeps for TP bypass %, as well as 2nd and 3rd regen stages


while True:
    regen = H2(h=turbo.h-dh_turbo, s=turbo.s)           # state 3
    pumped = H2(P=regen.P, s=start.s)                   # state 2
    dh_pump = pumped.h - start.h                        # empirical pump work
    q_nozzle = regen.h - pumped.h                       # required regenerative cooling
    err = dh_pump + dh_turbo                            # difference of turbine work and pump work
    print("Converging ---- enthalpy error:", err)
    if abs(err) < 1:
        print("Converged!")
        print("Turbopump enthalpy change:", dh_pump)
        print("Regen nozzle cooling:", q_nozzle)
        break
    dh_turbo -= err

bypass = H2(h=regen.h+q_regen1, P=regen.P)              # state 6
throt = H2(h=bypass.h, P=mix.P)                        # state 5



states = [start, pumped, regen, turbo, mix, bypass, throt, hot, cfe, pm, core]
flow1 = [start, pumped, regen, bypass, throt, mix, hot, cfe, pm, core]
flow2 = [regen, turbo, mix]
n = 5
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


plt.plot([point.s/1e3 for point in f1], [point.T for point in f1])
plt.plot([point.s/1e3 for point in flow2], [point.T for point in flow2])

T=range(700, 1200)
print("Max Pressure:", pumped.P)
s = H2.s.table.func_t(P1=pumped.P)
#plt.plot(s(T), T)
styles = '.*^v+x'
for i, point in enumerate(states):
    plt.plot(point.s/1e3, point.T, styles[i%6], label=i+1)


for point in states:
    print(point.T, "|", point.s)
plt.legend()
plt.ylabel("Temperature (K)")
plt.xlabel("Entropy (kJ/kg K)")

#current_values = plt.gca().get_yticks()
#plt.gca().set_yticklabels(['{:,.0f}'.format(x) for x in current_values])
#current_values = plt.gca().get_xticks()
#plt.gca().set_xticklabels(['{:,.0f}'.format(x) for x in current_values])

plt.show()


