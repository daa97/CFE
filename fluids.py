import numpy as np
import csv
from scipy.interpolate import CubicSpline, PPoly
import matplotlib.pyplot as plt
import os

debug = False
extra = False
def dbp(*args):
    "Prints arguments if in debug mode"
    if debug:
        print(*args)

class Table:
    '''Fluid property table which contains some data indexed to pressure (columns) and temperature (rows). '''

    def __init__(self, property, P: np.array, T: np.array, data: np.array):
        self.P_axis = P                     # pressure values corresponding to each column
        self.T_axis = T                     # temperature values corresponding to each row

        self.psplines = []                  # spline functions of pressure at tabulated temperatures
        self.tsplines = []                  # spline functions of temp at tabulated pressures
        self.psplines_Taxis = []            # temperature values corresponding to each spline f(P)
        self.tsplines_Paxis = []            # pressure values corresponding to each spline f(T)
        self.data = data                    # fluid property values
        self.tol = 1e-5                     # tolerance for removing duplicate solutions


        self.P_range = (min(self.P_axis), max(self.P_axis))     # pressure range of the table
        self.T_range = (min(self.T_axis), max(self.T_axis))     # temperature range of the table
        numdata = self.data[~np.isnan(self.data)]               # filter out NaN values for max and min
        self.data_range = (np.min(numdata), np.max(numdata))    # range of property values in the table
        self.property = property                                # metadata associated with the physical property


        # create spline fits as a function of pressure at specified temperatures
        for i in range(len(self.T_axis)):
            row = self.data[i,:]                # index values at each temperature
            x = self.P_axis[~np.isnan(row)]     # filter out pressure values which are missing property data
            y = row[~np.isnan(row)]             # filter out missing property data

            if len(y) > 1:                      # check if at least two data values exist in row
                # add spline function of pressure, and corresponding temperature value
                self.psplines.append(CubicSpline(x, y, extrapolate=extra))  
                self.psplines_Taxis.append(self.T_axis[i])
            
        # create spline fits as a function of temperature at specified pressures    
        for j in range(len(self.P_axis)):
            column = self.data[:,j]             # index values at each pressure
            x = self.T_axis[~np.isnan(column)]  # filter out temperature values which are missing property data
            y = column[~np.isnan(column)]       # filter out missing property data

            # check if at least two data values exist in column
            if len(y) > 1:
                # add spline function of temperature, and corresponding pressure value
                self.tsplines.append(CubicSpline(x, y, extrapolate=extra))
                self.tsplines_Paxis.append(self.P_axis[j])

    def linTP(self, T1,P1):
        row = np.nonzero(self.T_axis<=T1)[0][-1]
        col = np.nonzero(self.P_axis<=P1)[0][-1]
        Tlim = [self.T_axis[row], self.T_axis[row+1]]
        Plim = [self.P_axis[col], self.P_axis[col+1]]
        v0 = [[self.data[row,col], self.data[row,col+1]],
                        [self.data[row+1,col], self.data[row+1,col+1]]]
        v = np.array(v0)
        rowX = (T1 - Tlim[0])/(Tlim[1]-Tlim[0])
        colX = (P1 - Plim[0])/(Plim[1]-Plim[0])

        row_diff = (v[1,0]-v[0,0])*(1-colX) + (v[1,1]-v[0,1])*colX
        row_base = v[0,0]*(1-colX) + v[0,1] * colX
        # print(T1,Tlim,P1,Plim)
        # print(rowX, colX, row_diff, row_base)
        # print("*"*20)
        return row_diff * rowX + row_base

    def linT(self, Pval, Yval):
        col = np.nonzero(self.P_axis<=Pval)[0][-1]
        Plim = [self.P_axis[col+1], self.P_axis[col]]
        colX = (Pval - Plim[0])/(Plim[1]-Plim[0])
        Ypts = self.data[:,col]*(1-colX) + self.data[:,col]*colX
        for i in range(len(Ypts)-1):
            if min(Ypts[i], Ypts[i+1])<=Yval<=max(Ypts[i], Ypts[i+1]):
                break
        else:
            raise ValueError("Could not find temperature!")
        YX = (Yval - Ypts[i])/(Ypts[i+1] - Ypts[i])
        return self.T_axis[i] * (1-YX) + self.T_axis[i+1]*YX

    def check_range(self, P1=None, T1=None, Y1=None):
        "check if input property values are within range of table"
        self.check_single(P1, self.P_range, "Pressure")             # check pressure if provided
        self.check_single(T1, self.T_range, "Temperature")          # check temperature if provided
        self.check_single(Y1, self.data_range, self.property.name)  # check property value if provided

    def check_single(self, val, rang, name):
        "check if a provided value is finite and within provided range"
        if val is not None:
            if val<min(rang) or val>max(rang) or np.isnan(val):
                 raise ValueError(f"{name} value of {val} is not within tabulated range of {rang}")
    
    def interp(self, P1, T1) -> float:
        "Obtain property value from pressure and temperature through interpolation"
        self.check_range(P1=P1, T1=T1)      # ensure provided values are in table range
        FP = self.func_p(T1)                # create property function F(P) at given T
        FT = self.func_t(P1)                # create property function F(T) at given P
        Ans = (FP(P1) + FT(T1)) / 2         # average values to improve accuracy
        if np.isnan(Ans):                   # print warning if value is NaN
            print(f"WARNING: property '{self.property.name}' data not found at temperature {T1} & pressure {P1}")
        return Ans
    
    def func_t(self, P1) -> PPoly:
        "Property as a function of temperature at specified pressure P1"
        self.check_range(P1=P1)             # ensure provided values are in table range
        return self.func_1D(P1, self.psplines, self.psplines_Taxis)
    
    def func_p(self, T1) -> PPoly:
        "Property as a function of pressure at specified temperature T1"
        self.check_range(T1=T1)             # ensure provided values are in table range
        return self.func_1D(T1, self.tsplines, self.tsplines_Paxis)
    
    def func_1D(self, Z1, zsplines, x_axis) -> PPoly:
        '''Returns function y=F(x) at a specified z value.
        zsplines is a list of functions y(z) at x-intervals specified by x_axis'''

        y = np.array([func(Z1) for func in zsplines])      # get y(z) values for each x interval
        x = np.array(x_axis)[~np.isnan(y)]                    # filter out x-values where y is NaN
        y = y[~np.isnan(y)]                         # filter out NaN y-values
        return CubicSpline(x, y, extrapolate=False) # return interpolated spline function

    def get_p(self, T1, Y1) -> float:
        '''get pressure value from given temperature and property value (not recommended for enthalpy or cp tables).'''
        self.check_range(T1=T1, Y1=Y1)          # ensure provided values are in table range
        func = self.func_p(T1)                  # define 1D property function at given temperature
        return self.combine(func.solve(Y1))     # solve for pressure given property value
    
    def t_func_p(self,Y1) -> PPoly:
        '''Get temperature as a function of pressure for given property value'''
        self.check_range(Y1=Y1)             # ensure provided values are in table range        
        # solve for temperature at each pressure 
        T_points = []; P_points = []
        for i in range(len(self.tsplines)):
            tpt = self.tsplines[i].solve(Y1)
            if len(tpt)>0:
                T_points.append(self.combine(tpt))
                P_points.append(self.tsplines_Paxis[i])
        #T_points = [self.combine(F.solve(Y1)) for F in self.tsplines]
        # return function T(P) created from interpolation of these points
        return CubicSpline(P_points, T_points, extrapolate=False)
    
    def get_t(self, P1, Y1) -> float:
        '''Get temperature from pressure and property value'''
        self.check_range(P1=P1, Y1=Y1)          # ensure provided values are in table range
        func = self.func_t(P1)                  # define 1D property function at given pressure
        return self.combine(func.solve(Y1))     # solve for temperature given property value
        
    def combine(self, values: np.array) -> list:
        '''Remove duplicate and NaN solutions. 
        Should converge to a single solution or raise an error.'''
        vals = values[~np.isnan(values)]                # filter out NaN values
        unique = []                                 # list of unique solutions
        equal = lambda a,b: abs(a-b) < self.tol     # checks if two values are nearly equal to within specified tolerance
        for i, vi in enumerate(vals):               # loop through each value in list
            for j, vj in enumerate(vals[i:]):       # loop through values again to check pairs of values                
                if j!=0 and equal(vi, vj):          
                    break                           # don't include value i if equal to a later value
            else:
                unique.append(vi)      # keep value for output if it is unique
        if len(unique)!=1:                  # check that exactly 1 solution exists
            raise ValueError(f"No Unique Solution Found. Solutions list: {unique}")
        return unique[0]                    # return only solution

class FluidProperty:
    '''Representation of a fluid property such as pressure or entropy'''
    def __init__(self, tag:str, name:str, alt:list=[], rev:bool=True, file:str=None, units:str="", axis=False):
        self.name = name                                # full name of property
        self.isaxis = axis                              # specifies if property is pressure/temperature or another property
        self.tag = tag                                  # short tag to specify name
        self.units = units                              # units to display when showing the property
        self.file = file                                # file name of csv which contains property table
        self.reversible = rev                           # indicates whether the property can be used to determine state of the fluid
        self.known = (axis or (self.file is not None))  # indicates whether any data exists for this property
        # callnames lists all names which can be used to access the property
        self.callnames = [self.reformat(self.tag), self.reformat(self.name)] + [self.reformat(s) for s in alt]
        
    def read_table(self, fname=None):
        """reads property data table from file and assigns it to self.table"""
        if fname is not None:
            self.file = fname                       # update file name if new one is given
        if self.file is not None:
            self.known = True                       # indicate that there is an associated dataset for the property
            with open(self.file, mode='r', newline='') as f:    # read csv file, removing any thousands separators
                reader = csv.reader(f)
                data = np.array([[s.replace(',', '') for s in r] for r in reader])
            self.table = Table(self, P=data[0,1:].astype(float), T=data[1:,0].astype(float), data=data[1:,1:].astype(float))

    def knownas(self, string):
        "Tell whether a specific name can be used to access the property"
        return (self.reformat(string) in self.callnames) # return True if string is found in callnames

    def reformat(self, string):
        '''Formats string to avoid case sensitivity'''
        return string.upper().strip().replace(" ", "_")
    def __getstate__(self): 
        return self.__dict__
    def __setstate__(self, d):
        self.__dict__.update(d)

    def __getattr__(self, __name: str):
        '''allows direct access of table attributes from property'''
        if self.known and not self.isaxis:
            return self.table.__getattribute__(__name)
        else:
            raise AttributeError(f"{self}.{__name} not found")
            
    def __str__(self):
        """Controls printout behavior of property"""
        return f"[{self.tag}] {self.name} ({self.units})"
    
class FluidState:
    '''Representation of one particular state of a fluid with a given pressure, temperature, etc.'''
    def __init__(self, fluid, props, linear=False):
        self.properties = dict()
        self.fluid = fluid # Corresponding fluid with property tables
        self.linear = linear
        self.solve_properties(props) # solve for fluid properties
    def __getstate__(self): 
        return self.__dict__
    def __setstate__(self, d):
        self.__dict__.update(d)

    def __getattr__(self, __name: str):
        '''overrides default FluidState.attribute behavior for unrecognized attribute names. Allows parsing of any fluid property names desired'''
        for p in self.fluid.properties:
            if p.knownas(__name) and p.known:
                return self.properties[p]
    
    def solve_properties(self, props: dict):
        '''Solve for all properties of the fluid'''
        keys = list(props.keys())
        for prop in keys:
            if not prop.known:
                raise ValueError(f"Property '{prop}' values do not exist for this fluid")
            elif (not prop.reversible) and (not prop.isaxis):
                raise ValueError(f"Property '{prop}' cannot be used to solve for other properties")
        

        # Solve for both pressure and temperature if needed
        if not any([prop.isaxis for prop in keys]):
            self.double_solve(props)
        if self.fluid.T in keys:
            self.properties[self.fluid.T] = props[self.fluid.T]
        if self.fluid.P in keys:
            self.properties[self.fluid.P] = props[self.fluid.P]
        else:
            # solve for pressure
            other_prop = keys[keys is not self.fluid.T]
            self.properties[self.fluid.P] = other_prop.table.get_p( T1=self.T, Y1=props[other_prop])

        if not self.fluid.T in keys:
            # solve for temperature
            other_prop = keys[keys is not self.fluid.P]
            self.properties[self.fluid.T] = other_prop.table.get_t( P1=self.P, Y1=props[other_prop])
        
        # set value for each property
        for pr in self.fluid.properties:
            if pr in keys:
                self.properties[pr] = props[pr]
            elif pr.known and (not pr.isaxis):
                try:
                    if self.linear:
                        self.properties[pr] = pr.table.linTP(P1=self.P, T1=self.T)
                    else:
                        self.properties[pr] = pr.table.interp(P1=self.P, T1=self.T)
                except KeyError as err:
                    print("FAILED:", err, pr.name, "P,T:", self.P, self.T)
                    self.properties[pr] = 0
    def linsolve(props):
        keys = list(props.keys())
        vals = list(props.values())
        np.nonzero(keys[0].data>vals[0])

    def double_solve(self, props) -> "float, float":
        '''solve for both pressure and temperature as a function of other variables
        Should mostly only be called internally by other fluid state methods'''
        keys = list(props.keys())
        vals = list(props.values())
        t1_p = keys[0].t_func_p(vals[0])        # function T(P) given property 1
        t2_p = keys[1].t_func_p(vals[1])        # function T(P) given property 2
        
        pmax = min(max(t1_p.x), max(t2_p.x))    # upper bound of smaller range
        pmin = max(min(t1_p.x), min(t2_p.x))    # lower bound of smaller range
        x_new = []; c_new = []                  
        
        # iterate over pressure values for each
        for i in range(len(t1_p.x)):
            for j in range(len(t2_p.x)):
                p1 = t1_p.x[i]
                p2 = t2_p.x[j]
                if min(p1, p2)==pmax:
                    dbp("Last loop reached for property value intersection")
                    x_new.append(pmax)
                    break
                elif (max(p1,p2)>pmax) or (min(p1,p2)<pmin):
                    pass
                else:
                    # Append final x-value
                    if p1==p2:
                        c1 = t1_p.c[:,i]
                        c2 = t2_p.c[:,j]
                        x_new.append(p1)
                        c_new.append(c1 - c2)                
            else:
                continue
            break
        c=np.array(c_new).transpose(); x=np.array(x_new)
        difference = PPoly(c=c, x=x, extrapolate=False)
        if np.isnan(difference.solve(0)) or len(difference.solve(0))==0:
            a=1
        pressure_val = keys[0].combine(difference.solve(0))
        temperature = t1_p(pressure_val)
        self.properties[self.fluid.P] = pressure_val
        self.properties[self.fluid.T] = temperature

    # control printout behavior of fluid state
    def __str__(self):

        pstr = ["Fluid state with properties:"] + [f"{p}: {self.properties[p]}" for p in self.properties.keys()]
        return "\n\t".join(pstr)

class Fluid:

    '''Representation of a fluid medium. Contains data for the fluid across different states of varying pressure, temperature, etc.'''
    def __init__(self, name: str, prop_tables:'{str:str}'):
        self.name = name
        self.properties = [FluidProperty(tag='p', name='Pressure', alt=['PRESS'], rev=False, axis=True, units="Pa"),# Pressure
            FluidProperty(tag='t', name='Temperature', alt=['TEMP'], rev=False, axis=True, units="K"),              # Temperature
            FluidProperty(tag='h', name='Enthalpy', alt=['ENTH'], rev=True, units="J/kg"),                          # Enthalpy
            FluidProperty(tag='u', name='Internal Energy', rev=True, units="J/kg"),                                 # Internal Energy
            FluidProperty(tag='v', name='Specific Volume', alt=['VOLUME', 'VOL'], rev=True, units="m^3/kg"),        # Specific volume
            FluidProperty(tag='rho', name='Density', alt=['D'], rev=True, units="kg/m^3"),                          # Density
            FluidProperty(tag='s', name='Entropy', alt=['S', 'ENTROPY'], rev=True, units="J/kg K"),                 # Entropy
            FluidProperty(tag='cp', name='Specific Heat', alt=['C'], rev=False, units="J/kg K"),                    # Specific heat
            FluidProperty(tag='gamma', name='Specific Heat Ratio', alt=['Y', "GAMMAS"], rev=False),                 # gamma, sp. heat ratio
            FluidProperty(tag='mu', name='Dynamic Viscosity', alt=['VISC', 'VISCOSITY'], rev=False, units="Pa*s"),  # viscosity (dynamic)
            FluidProperty(tag='a', name='Speed of Sound',
                          alt=["SON_VEL", "SON_VELOCITY", "SONIC_VELOCITY"], rev=False, units="m/s"),               # speed of sound                     
            FluidProperty(tag='gibbs', name='Gibbs Free Energy', alt=['GFE'], rev=False, units="J/kg"),             # gibbs free energy
            FluidProperty(tag='k', name='Thermal Conductivity', alt=["CONDUCTIVITY"], rev=False, units="W/m K"),    # thermal conductivity
            FluidProperty(tag='m', name='Molecular Mass', alt=["MOLAR_MASS"], rev=False, units="kg/mol")]           # molecular mass
        for prop, fname in prop_tables.items():
            for p in self.properties:
                if p.knownas(prop):
                    p.read_table(fname)
    def process(self, state1, state2, linprops, n):
        '''Generates a list of `n` intermediate states given a start and ending state and the 
        names of two properties which should vary smoothly or stay constant. '''
        path = [state1]
        start = {linprops[0]:state1.__getattr__(linprops[0]), 
                 linprops[1]:state1.__getattr__(linprops[1])}
        final = {linprops[0]:state2.__getattr__(linprops[0]), 
                 linprops[1]:state2.__getattr__(linprops[1])}
        itp = lambda a,b,x: a*(1-x) + b*x       # function to interpolate properties
        for j in range(1,n):
            mid = dict()
            for k,v in start.items():
                mid[k] = itp(v,final[k],j/n)
            path.append(self.state(**mid))
        path.append(state2)
        return path
    
    def __str__(self):
        return f"Fluid object '{self.name}'"
    
    def __call__(self, **kwargs: float) -> FluidState:
        return self.state(**kwargs)
    def __getstate__(self): 
        return self.__dict__
    def __setstate__(self, d):
        self.__dict__.update(d)
    def __getattr__(self, __name: str):
        for p in self.properties:
            if p.knownas(__name):
                return p
        else:
            raise AttributeError(f"'Fluid' object has no attribute '{__name}'")
    
    def state(self, linear=False, **kwargs: float) -> FluidState:
        '''Returns a fluid state given two fluid properties. Property values can then be indexed based upon the given state.
        Examples:
        `pump_exit = H2.state(P=5, T=800); print(pump_exit.h)`
        `plt.plot(range(1, 5), [H2.state(pressure=Pi, s=1800).v for Pi in range(1, 5)])`
        '''
        if len(kwargs)==2:
            kw1 = list(kwargs.keys())[0]
            dbp("KWARGS:", kwargs)
            kw2 = list(kwargs.keys())[1]
        else:
            raise ValueError("Two properties are required to determine the fluid state")
        input_props = dict()

        for p in self.properties:
            if p.knownas(kw1):
                input_props[p] = float(kwargs[kw1])
            elif p.knownas(kw2):
                input_props[p] = float(kwargs[kw2])
        dbp("Input state properties:", input_props)
        if len(input_props)!=2:
            raise ValueError("Input properties not recognized")
        return FluidState(self, input_props, linear=linear)

class FluidsListClass:
    def __init__(self):
        maindir = "fluidprops" # Directory of csv property files
        self.prop_files = dict()
        for sub in [f.path for f in os.scandir(maindir) if f.is_dir()]:
            f = sub.split("\\")[-1]
            self.prop_files[f] = dict()
            for file in os.listdir(sub):
                tag = file.split(".")[0]
                self.prop_files[f][tag] = os.path.join(sub, file)
    def __getattr__(self, __name: str) -> Fluid:
        if __name=="H2":
            return Fluid("Hydrogen", self.prop_files["Hydrogen"])
        elif __name=="Air":
            return Fluid("Air", self.prop_files["Air"])


FluidsList = FluidsListClass()

