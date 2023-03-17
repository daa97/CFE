import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.interpolate as inter
from cfe_model import CFE
import matplotlib as mpl
import cfe_model as cm
from fluids import FluidsList

mpl.rc('font', family='Times New Roman',size="10")
mpl.rc('figure', figsize=(4.8,3.6))
mpl.rc('legend', labelspacing=0.05)
mpl.rc('savefig', dpi=800)
mpl.rc('lines', linewidth=1.2)
mpl.rc('axes', grid=True)
mpl.rc('grid', linewidth=0.25)
mpl.rc('mathtext', fontset="dejavuserif")
mpl.rc('xtick.minor', visible=True, size=1.5, width=0.5)
mpl.rc('ytick.minor', visible=True, size=1.5, width=0.5)
plt.rcParams['figure.constrained_layout.use'] =  True
np.set_printoptions(precision=5,suppress = True)

class labyrinth_seal:
    def __init__(self,labby_inputs):
        self.fluid = labby_inputs["fluid"]
        self.P_in = labby_inputs["Inlet pressure"]   
        self.P_out = labby_inputs["Outlet pressure"] 
        self.T_in = labby_inputs["Inlet temperature"]
        self.T_out = labby_inputs["Outlet temperature"]
        self.c = labby_inputs["Radial clearance"]
        self.w = labby_inputs["Tooth width"]
        self.s = labby_inputs["Tooth pitch"]
        self.h = labby_inputs["Tooth height"]
        self.R = labby_inputs["Shaft radius"] + self.h
        self.N = labby_inputs["Number of teeth"]
        self.A = 2 * np.pi * self.R * self.c
        self.check_lims()
        self.results = self.calculate_mass_flow()
        self.mdot = self.results[0]
        self.P_profile = self.results[1]
        self.Y_profile = self.results[2]
        self.rho_profile = self.results[3]

    def calculate_mass_flow(self):  
        i=0
        P = np.ones((self.N+1))
        P[0] = self.P_in
        P[-1] = self.P_out
        rho = np.ones((self.N+1))
        rho[1] = self.fluid(p = self.P_in, T = self.T_in).rho
        # print(rho[0])
        mdot = 0.1 * self.A * np.sqrt((self.P_in - self.P_out) * rho[1])
        mu = self.fluid(p=self.P_in,T=self.T_in).mu
        # print(mu)
        err = 1
        Y = np.ones((self.N))
        # print(P)
        # print(rho)
        while err>0.00001:    
            Re = mdot/ (np.pi * 2 * self.R * mu)
            g1 = (1-6.5*(self.c/self.s))
            g2 = 8.638*(self.c/self.s) *(self.w/self.s)
            g3 = 2.454*(self.c/self.s)
            g4 = 2.258*(self.c/self.s)
            ws = self.w/self.s
            gamma =(g1 - g2)*(Re+((g1)- g2)**(-1/((g3)+g4*ws**1.673)))**(g3+g4*ws**1.673)
            cd1 = (1.097-0.0029*self.w/self.c)*(1+(44.86*self.w/self.c)/Re)**(-0.2157)/np.sqrt(2)
            tdc = cd1*0.925*(gamma**0.861)
            P[1] = P[0] - (mdot/self.A)**2 / (2 * Y[0]**2 * tdc * cd1**2 *rho[1])
            Y[0] = 0.558 + 0.442*(P[1]/P[0])
            rho[2] = self.fluid(p=P[1], T=self.T_in).rho

            # print()
            # print("Reynolds number:",Re)
            # print("g1",g1)
            # print("g2",g2)
            # print("g3",g3)
            # print("g4",g4)
            # print("gamma:",gamma)
            # print("cd1:",cd1)
            # print("tdc:",tdc)
            for j in range(self.N-1):
                if j < 1:
                    pass
                else:
                    # print(j)
                    P[j+1] = P[j] - (mdot/self.A)**2 / (2 * Y[j]**2 * tdc**2 * rho[j+1])
                    # print("rho[j-1]:",rho[j-1])
                    # print("P[j+1]:",P[j+1])
                    rho[j+2] = self.fluid(p = P[j+1],t=self.T_in).rho
                    Y[j] = 0.558 + 0.442**(P[j+1]/P[j])
                    # print(P)
                    # print(rho)
            Y[-1] = 0.558 + 0.442**(P[self.N]/P[self.N-1])

            mdot_new = self.A*Y[self.N-1]*tdc*np.sqrt(2 * (rho[self.N]*(P[self.N-1]-P[self.N])))
            err = np.abs((mdot_new-mdot)/mdot)
            mdot = (mdot_new-mdot) * 0.05 + mdot
            # print("Mass flow:",mdot)
            # print("Error")
            i+=1
            if i > 2000:
                err = 0
                print("Could not converge")
        
        return [mdot,P,Y,rho]
    
    def check_lims(self):
        c_s = self.c/self.s
        print("c/s:",c_s)
        w_s = self.w/self.s
        print("w/s:",w_s)
        w_c = self.w/self.c
        print("w/c:",w_c)
        h_s = self.h/self.c
        print("h/s:",h_s)

    def plot_labby_geom(self):
        labby_x0 = 0
        labby_y0 = self.h + self.c
        labby_xs = np.zeros((self.N,))
        labby_ys = np.zeros((self.N,))
        labby_xs[0] = labby_x0
        labby_ys[0] = labby_y0

        for i in range(self.N):
            labby_xi = labby_x0 + i * self.s
            labby_yi = labby_y0
            labby_x = labby_xi
            labby_y = self.c
            plt.plot([labby_xi,labby_x],[labby_yi,labby_y],color="tab:blue")
            labby_xi = labby_x
            labby_yi = labby_y
            labby_x = labby_x + self.w
            plt.plot([labby_xi,labby_x],[labby_yi,labby_y],color="tab:blue")
            labby_xi = labby_x
            labby_yi = labby_y
            labby_y = labby_y + self.h
            plt.plot([labby_xi,labby_x],[labby_yi,labby_y],color="tab:blue")
            labby_xi = labby_x
            labby_yi = labby_y
            labby_x = labby_x + self.s - self.w
            plt.plot([labby_xi,labby_x],[labby_yi,labby_y],color="tab:blue")
        plt.plot([labby_x0,labby_x],[0,0],color = "r")
        plt.axis("scaled")
                  
        plt.show()     

def plot_labbys(teeth,press,PRs,air,p_given,plot_labby_inputs):      
        for i,num_teeth in enumerate(teeth):
            mdot_data = np.zeros((len(PRs),))
            for q,PR in enumerate(PRs):
                plot_labby_inputs["Number of teeth"] = num_teeth
                if p_given == "outlet":
                    plot_labby_inputs["Outlet pressure"] = press
                    plot_labby_inputs["Inlet pressure"] = PR * press
                
                elif p_given == "inlet":
                    plot_labby_inputs["Inlet pressure"] = press
                    plot_labby_inputs["Outlet pressure"] = press/PR
                print()
                print("Number of teeth:",num_teeth)
                labby = labyrinth_seal(plot_labby_inputs)
                print("Pressure ratio:",PR)
                mdot_data[q] = labby.mdot
                print("Mass flow rate:",mdot_data[q])
            plt.plot(PRs,mdot_data,label=f'{num_teeth} teeth')
        plt.legend()
        plt.title(f'Labyrinth seal mass flow for {press/6895:.2f} psia {p_given} pressure')
        plt.xlabel("Pressure ratio [-]")
        plt.ylabel("Mass flow rate [kg*s^-1]")
        plt.show()
                    
                
class stepped_labyrinth(labyrinth_seal):
    def __init__(self) -> None:
        pass

if __name__ == "__main__":
    air = FluidsList.Air
    labby_inputs = {
    "fluid" : air,
    "Radial clearance" : 0.001,# [m]
    "Tooth pitch" : 0.004, # [m]
    "Shaft radius" : 0.0471,
    "Tooth width" : 0.002,
    "Tooth height" : 0.004,
    "Number of teeth" : 10,
    "Inlet pressure" : 2*101325,
    "Outlet pressure" : 101325,
    "Inlet temperature" : 287.15,
    "Outlet temperature" : 287.15
}
