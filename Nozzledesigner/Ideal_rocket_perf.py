import numpy as np
from matplotlib import pyplot as plt
import math
from scipy.optimize import fsolve

class AonAStar:
    def value(self, M, gamma): # value of AonAstar for a given Mach and gamma
        return (1 / M) * ((1 + ((gamma - 1) / 2) * (M ** 2)) / ((gamma + 1) / 2)) ** ((gamma + 1) / (2 * (gamma - 1)))

    def mach(self, Arat, gamma): # mach number for a given area ratio
        def A_over_A_star(M):
            return (1 / M) * ((1 + ((gamma - 1) / 2) * (M ** 2)) / ((gamma + 1) / 2)) ** ((gamma + 1) / (2 * (gamma - 1))) - Arat

        subsonic_M = fsolve(A_over_A_star, 0.001)
        supersonic_M = fsolve(A_over_A_star, 5)
        #print('sub: ', subsonic_M, '\n super: ', supersonic_M)
        return subsonic_M, supersonic_M

def compflowtool(gamma,mach):
    p_0_over_p = (1+ ((gamma-1)/2)*mach**2)**(gamma/(gamma-1))
    T_0_over_T = (1+ ((gamma-1)/2)*mach**2)
    rho_0_over_rho = (1+ ((gamma-1)/2)*mach**2)**(1/(gamma-1))
    A_over_Astar = (1/mach)*((1+((gamma-1)/2)*mach**2)/((gamma+1)/2))**((gamma+1)/(2*(gamma-1)))
    MFP = (math.sqrt(gamma)*mach)/(1+((gamma-1)/2)*mach**2)**((gamma+1)/(2*(gamma-1)))

    return p_0_over_p,T_0_over_T,rho_0_over_rho,A_over_Astar,MFP

def rocket_perf(gamma, p_0, p_amb, T_0, g, epsilon, R, exit_mach):

    if exit_mach == []:
        exit_mach= AonAStar().mach(epsilon, gamma)[1]
    if epsilon ==[]:
        epsilon = AonAStar().value(exit_mach, gamma)
    p_0_over_p,T_0_over_T,rho_0_over_rho,A_over_Astar,MFP = compflowtool(gamma, exit_mach)
    p_e = p_0_over_p ** -1 * p_0
    T_e = T_0_over_T ** -1 * T_0

    c_F_idl = gamma * math.sqrt((2 / (gamma - 1)) * (2 / (gamma + 1)) ** ((gamma + 1) / (gamma - 1)) * (1 - (p_e / p_0) ** ((gamma - 1) / (gamma)))) + (p_e / p_0 - p_amb / p_0) * epsilon
    c_star_idl = math.sqrt(1 / gamma * ((gamma + 1) / 2) ** ((gamma + 1) / (gamma - 1)) * R * T_0)

    mdotoverA = MFP * p_0 / (math.sqrt(R * T_0))
    u_e = exit_mach * math.sqrt(gamma * R * T_e)
    u_eq = u_e + (p_e - p_amb) / mdotoverA
    I_sp = u_eq / g
    Tau_over_A_t = c_F_idl * p_0
    return c_F_idl, c_star_idl, I_sp, Tau_over_A_t, u_eq, u_e
