import numpy as np
import matplotlib as plt
from testforhs.py import test_turb

if __name__ == "__main__":
    self = test_turb
    def make_hub_and_shroud():
        n = 2
        N_p = 100 # Number of points used to make the hub and shroud curves
        if self.r_4-self.r_h5 > self.z_r:
            R_c = self.z_r
            L = self.r_4 - self.z_r
        else:
            R_c = self.r_4-self.r_h5
            L = self.z_r - (self.r_4 - self.r_h5)
        hub_points = np.zeros(N_p,2)
        np.pi/2 / N_p
        for i,theta in enumerate(range(0, np.pi/2 / N_p, np.pi/2)):
            hub_points[i,0] = R_c * np.sin(theta)
            hub_points[i,2] = R_c * np.cos(theta)
        plt.plot(hub_points)
        plt.show()