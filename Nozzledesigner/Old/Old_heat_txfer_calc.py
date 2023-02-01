from Contour_designer import *
from Ideal_rocket_perf import *

#Constants
gamma = 1.4 # specific heat ratio |  This is not necessarily true and needs to be solved for with CEA!
MW = 2 # This is not necessarily true and needs to be solved for with CEA!
R = 8314/MW

# Overall Rocket Inputs

g = 9.81
F = 5e3 #N thrust
epsilon = 120 #expansion area ratio (Ae/At) or Ae/A*)

# "Chamber" Properties
p_c = 550*6894.76 #Pa chamber pressure (stagnation pressure before entering nozzle)
p_amb = 0 #ambient pressure (vacuum)
T_c = 3800 #K
rho_c = p_c/(R*T_c)


subsonic_M, M = a_on_a_star(gamma, epsilon)
C_f, c_star_idl, I_sp, Tau_over_A_t, u_eq, u_e, = rocket_perf(gamma, p_c, p_amb, T_c, g, epsilon, R, M)
A_t = F/C_f/p_c
d_t = math.sqrt(4*A_t/math.pi)
A_e = A_t*epsilon
d_e = math.sqrt(4*A_e/math.pi)

# Nozzle Geometry Starting Properties

L_star = .1 # characteristic combustion chamber length
alpha = math.radians(15) #half angle cone nozzle used to initially design this one; this is default
theta_n = math.radians(30) # starting angle for parabola approximation (see Huzel and Huang)
theta_e = math.radians(8) # parabola ending angle (see Huzel and Huang)
cone_noz_perc = .8 # percentage of equivalent cone nozzle used to initially design the bell contour (this is default)
data_pts = 100

# Form Nozzle Contour

X,Y = contour(L_star,alpha,epsilon,A_t, theta_n, theta_e, cone_noz_perc,data_pts,showplot=True)

#Create array of area ratios throughout nozzle
areas = np.pi*np.array(Y)**2
area_rats = areas/A_t

throat_idx = np.argmin(area_rats)
ic(throat_idx)
M_subsonic = []
M_supersonic = []

#Create array of mach numbers from area ratios
for i in area_rats[0][:]:
    M_sub = a_on_a_star(gamma, i)[0]
    np.append(M_subsonic, M_sub)

for i in area_rats[1][:]:
    M_super = a_on_a_star(gamma, i)[1]
    np.append(M_supersonic,M_super)

M_subsonic = np.append(M_subsonic,np.array(1))
M_dist = np.concatenate((M_subsonic, M_supersonic))

print(M_dist)
T_dist = []
V_dist = []
h_g_dist = []
M_dist = np.array(M_dist)

# create 1D thermo properties from mach number distribution
for i in M_dist:
    TO_over_T= compflowtool(gamma, i)[2]
    rho_O_over_rho = compflowtool(gamma, i)[3]

    T = T_c/TO_over_T
    V = np.sqrt(gamma*R*T)*i
    rho = rho_c/rho_O_over_rho
    h_g = (rho*V)**.8

    T_dist = np.append(T_dist, T)
    V_dist = np.append(V_dist, V)
    h_g_dist = np.append(h_g_dist,h_g)

T_wg = 3800 #K


# Solve for nozzle heat transfer
ic(np.shape(area_rats))
ic(np.shape(T_dist))
plt.plot(area_rats,T_dist)
plt.show()
plt.legend()


#using most basic approximation
q_dot_v = np.sum(h_g_dist*(np.array(T_dist)-np.array(T_wg)))
ic(q_dot_v)

#next-level approximation:

# Determine temperature given an area ratio

#(convert the following to for loops later)
temps_sub= rcea.sub.temps(area_rat,p_c,T_c) # temperature in nozzle calculated by CEA for a given subsonic area ratio
temp_th = rcea.throat.temp(area_rat,p_c,T_c)
temps_sup = rcea.sup.temps(area_rat,p_c,T_c)

T_dist_CEA = [temps_sub, temp_th, temps_sup]

# heat transfer calculation

T_L = 50 #K #Specify starting bulk coolant temperature
t_w= .1 #m wall thickness
k_w = 1 #thermal conductivity of wall material (copper?)  in SI units
h_L = 1 #coolant convective coefficients (hopefully find somewhere)

q_dot_v = (temps_sub - T_L)/(1/h_g_dist +t_w/k_w +1/h_L)




