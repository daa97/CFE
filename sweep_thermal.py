import prop_vary as pv
import thermal_profile as tp

vary = pv.pcore
#vary = pv.mdot
#vary = pv.N
#vary = pv.L

num_cells = 400
n_pts = 50

tp.sweep(vary, num_cells, n_pts)