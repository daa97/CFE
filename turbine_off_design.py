import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.interpolate as inter
import csv
import time
from fluids import FluidsList
H2 = FluidsList.H2
air = FluidsList.Air
import matplotlib as mpl
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

class off_design_analysis:
    def __init__(self,turbine,stator,op_point,i):
        self.DPT = turbine
        DPT_geom = self.DPT.rotor_geom
        self.r_4 = DPT_geom["inlet radius"]
        self.b_4 = DPT_geom["inlet blade height"]
        self.beta_4 = np.pi/2 - DPT_geom["inlet blade angle"]
        self.r_5h = DPT_geom["outlet hub radius"]
        self.r_5s = DPT_geom["outlet shroud radius"]
        self.b_5 = DPT_geom["outlet blade height"]
        self.beta_5 = np.pi/2 - DPT_geom["outlet blade angle"]
        self.beta_5s = DPT_geom["outlet shroud blade angle"]
        self.beta_5h = DPT_geom["outlet hub blade angle"]
        self.z_r = DPT_geom["axial length"]
        self.t_lead = DPT_geom["leading edge thickness"]
        self.t_lead = DPT_geom["leading edge thickness"]
        self.t_trail = DPT_geom["trailing edge thickness"]
        self.n_r = op_point["number of blades"]
        self.r_3 = DPT_geom["stator outlet radius"]

        self.DPS = stator
        DPS_geom = self.DPS.geom
        self.operating_point = op_point
        self.mass_guess = self.calc_mass_flow_ns
        
    # def calc_mass_flow_ns(self):
    #     self.DPT.
class component:
    def __init__(self,geometry_matrix, angle_convention,component_type) -> None:
        """Component type - Dict that tells what component it is, 
        i.e. rotor or stator, whether is is rotating or stationary, """
        self.radii = geometry_matrix[0]
        self.z = geometry_matrix[1]
        self.b = geometry_matrix[2]
        self.o = geometry_matrix[3]
        self.beta = geometry_matrix[4]
        self.theta = geometry_matrix[5]
        self.m = geometry_matrix[6]
        self.L = geometry_matrix[7]

    def change_angle_convention(self,angles):
        if self.angle_convention.lower() == "baines":
            old_angles = angles
            new_angles = (angle - np.pi/2 for angle in old_angles)




class station: 
    """Not sure I'll need this but if I do it's here. 
    Probably will be nice to assign 3 stations to each 
    component and make it easy to find values for BL losses"""
    def __init__(self) -> None:
        pass

class rotor(component):
    def __init__(self) -> None:
        pass

class nozzle(component):
    def __init__(self) -> None:
        pass

class vaneless_passage:
    def __init__(self) -> None:
        pass
def boundary_losses(component):
    pass

def performance_estimation_aungier(num_components,mass_flow_guess="calculated",p_dis = 101325):
    pass