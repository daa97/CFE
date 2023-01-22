import numpy as np
import matplotlib as plt
from testforhs import init_turb

if __name__ == "__main__":

    def make_hub_and_shroud():
        n = 2
        N_p = 100 # Number of points used to make the hub and shroud curves
        if init_turb.r_4-init_turb.r_h5 > init_turb.z_r:
            R_c = init_turb.z_r
            L = init_turb.r_4 - init_turb.z_r
        else:
            R_c = init_turb.r_4-init_turb.r_h5
            L = init_turb.z_r - (init_turb.r_4 - init_turb.r_h5)
        hub_points = np.zeros(N_p,2)
        np.pi/2 / N_p
        for i,theta in enumerate(range(0, np.pi/2 / N_p, np.pi/2)):
            hub_points[i,0] = R_c * np.sin(theta)
            hub_points[i,2] = R_c * np.cos(theta)
        plt.plot(hub_points)
        plt.show()

    init_turb.make_hub_and_shroud()