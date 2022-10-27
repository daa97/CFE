from fluids import *

import matplotlib.pyplot as plt

q_nozzle = 1e4
dP_pump = 10
dh_turbo = -4e3 
y_turb = 0.8            # fraction going into turbopump
q_regen1 = 2e3          # bypass heating states 3-6
q_regen2 = 5e3          # pre-cfe heating states 5-8
dh_cfeturb = -2e3
dp_pm = -.09
q_uran = 3e4

y_throt = 1-y_turb
H2 = Fluid("Hydrogen", prop_files)

start = H2(P=0.1, T=200)                                # state 1
pumped = H2(s=start.s, P=start.P+dP_pump)               # state 2
regen = H2(h=pumped.h+q_nozzle, P=pumped.P)             # state 3
turbo = H2(h=regen.h+dh_turbo, s=regen.s)               # state 4
bypass = H2(h=regen.h+q_regen1, P=regen.P)              # state 6
throt = H2(h=bypass.h, P=turbo.P)                        # state 7
mix = H2(h=y_turb*turbo.h+y_throt*throt.h, P=turbo.P)   # state 5
hot = H2(h=mix.h+q_regen2, P=mix.P)                     # state 8
cfe = H2(h=hot.h+dh_cfeturb, s=hot.s)                   # state 9
pm = H2(h=cfe.h, P=cfe.P+dp_pm)                         # state 10
core = H2(h=pm.h+q_uran, P=pm.P)                        # state 11

flow1 = [start, pumped, regen, bypass, throt, mix, hot, cfe, pm, core]
flow2 = [start, pumped, regen, turbo, mix, hot, cfe, pm, core]


states = [start, pumped, regen, turbo, bypass, throt, mix, hot, cfe, pm, core]


plt.plot([point.s for point in flow1], [point.T for point in flow1])
plt.plot([point.s for point in flow2], [point.T for point in flow2])
plt.plot([point.s for point in states], [point.T for point in states], '.')
plt.show()


