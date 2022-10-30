from fluids import *

import matplotlib.pyplot as plt
def moment(mass, R1, R2):
    return (1/2) * mass * (R1**2 + R2**2)

def mass(R1, R2, length, density):
    return np.abs(np.pi * (R2**2 - R1**2) * length) * density

l = 0.84
r1 = 0.03
r2 = 0.045
rho_u = 17100
void_frac = 0.5
omega = 7000 * (2 * np.pi / 60)
t_init = 10 # seconds for spin-up

rho_uh = (1 - void_frac) * rho_u
mass_u = mass(r1, r2, l, rho_uh)
moment_u = moment(mass_u, r1, r2)
ke_cfe = (1/2) * moment_u * omega**2
cfe_power = ke_cfe / t_init
cfe_flow = 0.1

P_tank = 1.5e5
T_tank = 40
p_core = 5e6
T_core = 5000
T_pm = 800
dp_uran = -2.2e6
dp_pm = -9e3

eta_turbopump =1
eta_cfeturb =1


dh_cfeturb = - cfe_power / cfe_flow
q_uran = 3e7

dh_turbo = -6e6
q_nozzle = 8e5

y_turb = 1            # fraction going into turbopump
dh_pump = -dh_turbo * eta_turbopump * y_turb

q_regen1 = 1e6 *0         # bypass heating states 3-6
q_regen2 = 1e6 *0         # pre-cfe heating states 5-8
y_throt = 1-y_turb

H2 = Fluid("Hydrogen", prop_files)

core = H2(T=T_core, P=p_core)                           # state 11
pm = H2(T=T_pm, P=core.P - dp_uran)                     # state 10
cfe = H2(h=pm.h, P=pm.P - dp_pm)                        # state 9
hot_isent = H2(h=cfe.h-dh_cfeturb, s=cfe.s)             # state 8
hot = H2(h=cfe.h-dh_cfeturb*eta_cfeturb, P=hot_isent.P) 
mix = H2(h=hot.h-q_regen2, P=hot.P)                     # state 8
turbo = H2(h=mix.h, P=mix.P)                # state 4
start = H2(P=P_tank, T=T_tank)  
print(start)
#tfunc1 = H2.s.table.t_func_p(Y1=start.s)
#tfunc2 = H2.s.table.t_func_p(Y1=turbo.s)
#regen = H2(h=turbo.h-dh_turbo, s=turbo.s)               # state 3

while True:
    regen = H2(h=turbo.h-dh_turbo, s=turbo.s)               # state 3
    pumped = H2(P=regen.P, s=start.s)
    hp = pumped.h - start.h
    res = hp + dh_turbo
    print("RESIDUAL:", res)
    if abs(res) < 10:
        break
    dh_turbo -= res


# hp + dh_turbo = 0
# pumped.h  = start.h - dh_turbo
# 
 
#pumped = H2(s=start.s, h=start.h+dh_pump)                # state 2
#regen = H2(h=pumped.h+q_nozzle, P=pumped.P)             # state 3
#turbo = H2(h=regen.h+dh_turbo, s=regen.s)               # state 4
bypass = H2(h=regen.h+q_regen1, P=regen.P)              # state 6
throt = H2(h=bypass.h, P=mix.P)                        # state 5







states = [start, pumped, regen, turbo, bypass, throt, mix, hot, cfe, pm, core]
print(len(states))
flow1 = [start, pumped, regen, bypass, throt, mix, hot, cfe, pm, core]
flow2 = [start, pumped, regen, turbo, mix, hot, cfe, pm, core]
#plt.plot([point.s for point in flow1[:-1]], [point.T for point in flow1[:-1]])
plt.plot([point.s for point in flow2[:-1]], [point.T for point in flow2[:-1]])
i=11
T=range(700, 1200)
print("Max Pressure:", pumped.P)
s = H2.s.table.func_t(P1=pumped.P)
#plt.plot(s(T), T)
for point in states[:-1]:
    i -=1
    plt.plot(point.s, point.T, '.', label=i)


for point in states:
    print(point.T, "|", point.s)
plt.show()


