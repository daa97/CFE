import math
from matplotlib import pyplot as plt
import numpy as np
from icecream import ic

A_t =
L_star =
theta_n =
theta_e =
cone_noz_perc =
alpha =
data_pts =

p_c =
T_c =
e_c =


class Contour:
    def __init__(self, p_c, T_c, A_t, e, start_geometry):
        self.A_t = A_t
        self.p_c = p_c
        self.T_c = T_c
        self.e = e
        self.e_c = e_c
        self.start_geometry = start_geometry
        '''start_geometry = {"L_star": L_star, "Theta_n": theta_n, "alpha": alpha, "theta_e": theta_e,
                          "cone_noz_perc": cone_noz_perc, "data_pts": data_pts, "e_c": e_c}'''

    def make_contour(self,showplot=False):

        #derived from inputs
        d_t =math.sqrt(self.A_t*4/math.pi)
        r_t = d_t/2
        r_e = np.sqrt(epsilon)*r_t
        lmbda = (1+math.cos(alpha))/2

        def cone_length(alpha,epsilon,r_t,r_e):
            return (r_t * (math.sqrt(self.e) - 1) + r_e * (1/math.cos(self.alpha) - 1)) / math.tan(self.alpha)

        #Converging section
        V_c = self.start_geometry["L_star"]*self.A_t

         #(rationale: typical values are 2-5, this needs optimization)
        d_c = math.sqrt(self.e_c)*d_t
        ic(d_c)
        r_c = d_c/2
        conv_cone_angle = math.radians(20)
        conv_cone_len = cone_length(math.radians(20), self.e_c,r_t,r_c)
        conv_cone_vol = math.pi/3*conv_cone_len*(r_c**2+r_t**2+r_c*r_t)

        V_c_cyl = V_c-conv_cone_vol
        L_c_cyl = V_c_cyl/(self.e_c*A_t)

        #Diverging section
        L_n = cone_noz_perc* cone_length(math.radians(15), self.e, r_t, r_e)

        #N coordinates (point where it transitions to parabola
        N_t = .385*r_t*math.sin(math.radians(30))
        ic(N_t)

        N_a = r_t +.385*r_t*(1-math.cos(theta_n))
        E_t = L_n
        E_a = r_e

        # Region I

        r1 = 1.5*r_t
        th1 = np.linspace(-math.pi,-math.pi/2, int(math.floor(data_pts/3)))
        x1 = r1*np.cos(th1)
        y1 = r1*np.sin(th1)+2.5*r_t


        # Region II
        def n_point_finder(r_t, theta_n):
            r = .385 * r_t
            slope = math.tan(theta_n)
            theta = np.linspace(-math.pi/2, 0, int(math.floor(data_pts/3)))
            theta.astype(int)
            x = r * np.cos(theta)
            y = r * np.sin(theta)+1.385*r_t
            grad_array = np.gradient(y, x)

            def find_nearest(array, value):
                idx = np.searchsorted(array, value, side="left")
                if idx > 0 and (idx == len(array) or math.fabs(value - array[idx - 1]) < math.fabs(value - array[idx])):
                    return idx - 1
                else:
                    return idx

            idx = find_nearest(grad_array, slope)
            x_n = x[idx]
            y_n = y[idx]
            return x_n, y_n


        x_n, y_n = n_point_finder(r_t,theta_n)
        R = .385*r_t
        th2 = np.linspace(-np.pi/2,(theta_n-np.pi/2), int(math.floor(data_pts/3)))
        x2 = R*np.cos(th2)
        y2 = R*np.sin(th2)+1.385*r_t

        # Region III
        # fit parabola between two points
        y_e = r_e
        x_e = L_n
        m_n = math.tan(theta_n)
        m_e = math.tan(theta_e)

        A_mat = np.array([[2*y_n*m_n, m_n],
                     [2*y_e*m_e, m_e]])
        B_mat = np.array([1, 1])
        A,B = np.linalg.solve(A_mat, B_mat)
        C = x_n -A*y_n**2 -B*y_n

        y3 = np.array(np.linspace(y_n,y_e,int(math.floor(data_pts/3))))
        x3 = A*y3**2+B*y3+C

        # Plots
        plt.plot(x1, y1, label='1')
        plt.plot(x2, y2, label='2')
        plt.plot(x3, y3, label='3')
        plt.axis('square')
        if showplot:
            plt.show()
        X = [x1, x2, x3]
        Y = [y1,y2,y3]
        self.contour = X,Y

