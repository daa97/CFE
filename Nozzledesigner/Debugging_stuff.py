'''def eq2iteratethrough(T_wg, Nozzle, mu_hg_dist, cp_hg_dist, gamma_dist, M_dist, h_L_dist, k_w, t_w, i):
    ## GO BACK AND FIND OUT WHAT w should BE!!
    # This equation below comes from assuming 1-D steady-state heat transfer with negligible radiation heat transfer
    eq = bartz(d_t=Nozzle.d_t, mu_c=mu_hg_dist[i], Cp_c=cp_hg_dist[i], Pr_c=Pr_hg_dist[i], p_c=Nozzle.p_c,
               T_c=Nozzle.T_c,
               Cstar=cstar_dist[0], r_c=Nozzle.r_c, A=Nozzle.area_rats[i], w=.4, gamma=gamma_dist[i], M=M_dist[i],
               T_w=T_wg) * (T_hg_dist[i] - T_wg) + k_w / t_w * (
                     T_wg - ((k_w / t_w) * T_wg + h_L_dist[i] * T_L[i]) / (h_L_dist[i] - k_w / t_w))
    return eq


T_wg = scipy.optimize.brentq(
    eq2iteratethrough(i=i, Nozzle=Nozzle, mu_hg_dist=mu_hg_dist, cp_hg_dist=cp_hg_dist, gamma_dist=gamma_dist,
                      M_dist=M_dist, k_w=k_w, t_w=t_w, h_L_dist=h_L_dist), 0, 1e8)'''


import scipy.optimize
import icecream as ic

def eq2solve(x,a,b):
    eq = -x**2 +a*b
    return eq
a=3
b=3

solution = scipy.optimize.brentq(eq2solve,0,10,args=(a,b))
print(solution)