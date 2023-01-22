import math
from Ideal_rocket_perf import *


#constants
g = 9.81

# inputs
F = 5e3 #N thrust
p_c = 550*6894.76 #Pa chamber pressure (stagnation pressure before entering nozzle)
p_amb = 0 #ambient pressure (vacuum)
gamma = 1.4 # specific heat ratio
epsilon = 25 #expansion area ratio (Ae/At) or Ae/A*)
T_c = 5500 #K
Rbar = 8314 #j/kgK
MW = 2
L_star = 35*.0254 #m

#def nozzle_designer_idl(F, p_c, C_f, epsilon, L_star):
    # region I
subsonic_M, M = a_on_a_star(gamma, epsilon)
C_f, c_star_idl, I_sp, Tau_over_A_t, u_eq, u_e, = rocket_perf(gamma, p_c, p_amb, T_c, g, epsilon, Rbar/MW, M)
A_t = F/C_f/p_c
d_t = math.sqrt(4*A_t/math.pi)
A_e = A_t*epsilon
d_e = math.sqrt(4*A_e/math.pi)









#=========================================
# region II
r_t = d_t/2
L_c = .05*math.sqrt(d_t) # Conditional Length
epsilon_c = L_star/L_c # contraction area ratio
A_c = A_t*epsilon_c # Chamber area
d_c = math.sqrt(4*A_c/math.pi)
r_c = d_c/2

# nozzle contour using Vitoshinsky's formula (see reference 10)
'''x = np.linspace(0.00,4.00, num =100)
y = r_t / math.sqrt(
(1 - (1 - (r_t / r_c) ** 2) * (1 - (x / 1.5 / r_c) ** 2) ** 2) / (1 - ((x / 1.5 / r_c) ** 2) / 3) ** 3)'''

#Starter nozzle contour using sutton method


V_II = .013 # create function that finds volume of the above contour
V_cc = A_t * L_star
V_I = V_cc - V_II
L_I = V_I/A_c

#return A_t, d_t, A_e, d_e, x, y
print("d_t= ", d_t, "\n d_e= ", d_e, "\n A_t = ", A_t, "I_sp= ", I_sp)


'''def heat_transf_calc
    T_aw = T_c*
subsonic_M, M_e = a_on_a_star(gamma, epsilon)
R = Rbar/MW
c_F_idl, c_star_idl, I_sp, Tau_over_A_t, u_eq, u_e, I_sp = rocket_perf(gamma, p_c, p_amb, T_c, g, epsilon, R, M_e)
m_dot = F/I_sp/g
'''
