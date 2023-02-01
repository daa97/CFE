
from fluids import *
from icecream import ic
import math
#https://www.htflux.com/en/documentation/boundary-conditions/surface-resistance-heat-transfer-coefficient/heat-transfer-of-pipe-flows/
# https://www.engineeringtoolbox.com/conductive-heat-loss-cylinder-pipe-d_1487.html
mdot = 2.2
OD = 1.315 * .0254
ID = 1.049 * .0254
t = .133 * .0254  # schedule 40 pipe thickness
def pipe_heat_xfer(T,P,mdot,OD,ID,t,T_wh,k_w):
    nz_cool_start = H2(T=T,P=P)
    fluid_IC = nz_cool_start

    A_c = np.pi*(ID)**2/4 # assume 1" pipe diameter for now
    char_len = 4*ID #using 4X the pipe diameter as the characteristic length

    k_L = fluid_IC.k
    rho_L = fluid_IC.rho
    v_L = mdot / (rho_L * A_c)

    mu_L = fluid_IC.mu
    Re = rho_L * v_L * char_len / mu_L
    cp_L = fluid_IC.cp
    Pr = mu_L * cp_L / k_L

    h_L = k_L / char_len * .023 * Re ** .8 * Pr ** .4  # Dittus Boelter correlation for turblent flow pipes
    T_L = fluid_IC.t
    ic(T_L)

    q = 2*np.pi* (T_L -T_wh) / (math.log((OD/2)/(ID/2)) / k_w) # very basic correlation for heat transfer for pipe length
    return q
    ic(q)