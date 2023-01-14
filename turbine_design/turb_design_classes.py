import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.interpolate as inter
import csv
from fluids import *

class IdealFluid:
    """Defines a fluid whose properties can be calculated using the ideal gas law--used for model validation"""
    def __init__(self,cp,gamma):
        self.cp = cp
        self.gamma = gamma
        self.R = self.cp * (1 - 1/self.gamma)

class IdealFluidState:
    def __init__(self,idealFluid,**kwargs: float):
        self.fluid = idealFluid
        self.t = kwargs["t"]
        self.p = kwargs["p"]
        self.rho = self.find_rho()

    def find_rho(self):
        rho = self.p / (self.fluid.R * self.t)
        return rho

class RadialTurbine:
    def __init__(self,CFE,static_turb_inputs,dynamic_turb_inputs,i):
        self.i = i
        self.CFE = CFE
        self