from turb_design_classes import *

if __name__ == "__main__":
    static_turb_inputs_t100 = {
        "T_01" : 1056.5, #[K] Stagnation (total) inlet temperature
        "P_01" : 0.5804, #[MPa] Stagnation (total) inlet temperature #FIXME - change to MPa and outlet pressure for CFE
        "mass_flow" : 0.33,
        "PR_ts" : 5.73,
        "N_r" : 16,
        "N_n" : 19,
        "W_dot" : 121000,
        "C_m1" : 0,
        "C_theta1" : 0
    }