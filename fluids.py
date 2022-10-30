from re import A
import numpy as np
import csv
from scipy.interpolate import CubicSpline, PPoly
import matplotlib.pyplot as plt
import os.path

debug = False
extra = False
def dbp(*args):
    "Prints arguments if in debug mode"
    if debug:
        print(*args)

class Table:
    '''Fluid property table which contains some data indexed to pressure (columns) and temperature (rows). '''

    def __init__(self, property, P: np.array, T: np.array, data: np.array):
        self.P_axis = P
        self.T_axis = T
        self.csp_axis = []
        self.cst_axis = []
        self.data = data            # fluid property values
        self.tol = 1e-5             # tolerance for removing duplicate values
        self.P_range = (min(self.P_axis), max(self.P_axis))
        self.T_range = (min(self.T_axis), max(self.T_axis))
        self.data_range = (np.min(self.data), np.max(self.data))
        self.property = property


        # create spline fits as a function of pressure at specified temperatures
        self.cs_p = [] 
        for i in range(len(self.T_axis)):
            row = self.data[i,:]
            x = self.P_axis[~np.isnan(row)]
            y = row[~np.isnan(row)]
            if len(y) > 1:
                self.cs_p.append(CubicSpline(x, y, extrapolate=extra))
                self.csp_axis.append(self.T_axis[i])
            
        # create spline fits as a function of temperature at specified pressures    
        self.cs_t = []                   
        for j in range(len(self.P_axis)):
            column = self.data[:,j]
            x = self.T_axis[~np.isnan(column)]
            y = column[~np.isnan(column)]
            if len(y) > 1:
                self.cs_t.append(CubicSpline(x, y, extrapolate=extra))
                self.cst_axis.append(self.P_axis[j])
    def check_range(self, P1=None, T1=None, Y1=None):
        "check if input property values are within range of table"
        if P1 is not None:
            if P1<self.P_range[0] or P1>self.P_range[1]:
                raise ValueError(f"Pressure of {P1} is not within tabulated range of {self.P_range}")
        if T1 is not None:
            if T1<self.T_range[0] or T1>self.T_range[1]:
                raise ValueError(f"Temperature of {T1} is not within tabulated range of {self.T_range}")
        if Y1 is not None:
            if Y1<self.data_range[0] or Y1>self.data_range[1]:
                raise ValueError(f"Property value of {Y1} is not within tabulated range of {self.data_range}")

    def interp(self, P1, T1) -> float:
        "Obtain property value from pressure and temperature"
        self.check_range(P1=P1, T1=T1)
        FP = self.func_p(T1)
        FT = self.func_t(P1)
        Ans = (FP(P1) + FT(T1)) / 2
        if np.isnan(Ans):
            print(f"WARNING: fluid property '{self.property.name}' could not be found at temperature {T1} & pressure {P1}")
        return Ans
    
    def func_t(self, P1) -> PPoly:
        "Property as a function of temperature at specified pressure P1"
        self.check_range(P1=P1)
        at_p1 = []
        tvals = []
        for i in range(len(self.csp_axis)):
            yval = self.cs_p[i](P1)
            if np.isfinite(yval):
                at_p1.append(yval)
                tvals.append(self.csp_axis[i])

        return CubicSpline(tvals, at_p1, extrapolate=extra)
    
    
    def func_p(self, T1) -> PPoly:
        "Property as a function of pressure at specified temperature T1"
        self.check_range(T1=T1)
        at_t1 = []
        pvals = []
        for j in range(len(self.cst_axis)):
            yval = self.cs_t[j](T1)
            if np.isfinite(yval):
                at_t1.append(yval)
                pvals.append(self.cst_axis[j])


        return CubicSpline(pvals, at_t1, extrapolate=extra)
    
    
    def get_p(self, T1, Y1) -> float:
        '''get pressure value from given temperature and property value (not recommended for enthalpy or cp tables).'''
        self.check_range(T1=T1, Y1=Y1)
        func = self.func_p(T1)
        solns = self.combine(func.solve(Y1))
        if len(solns)!=1:
            print(solns)
            raise ValueError("No Unique Solution Found")
        return solns[0]
    
    def t_func_p(self,Y1) -> PPoly:
        '''Get temperature as a function of pressure for given property value'''
        self.check_range(Y1=Y1)
        # find temperature at tabulated pressures and given property value Y1
        intercepts = dict()
        dbp(Y1)
        dbp(self.T_axis)
        for i in range(len(self.cs_t)):
            Ti = self.combine(self.cs_t[i].solve(Y1))
            dbp(Ti)
            # Check that there is only one valid temperature solution at each pressure
            if len(Ti)==0:
                pass
            elif len(Ti)==1:
                intercepts[self.cst_axis[i]] = Ti[0]
            else:
                raise ValueError("No Unique Solution Found")
        
        # return function T(P) created from interpolation of these points
        dbp(list(intercepts.keys()), list(intercepts.values()))
        return CubicSpline(list(intercepts.keys()), list(intercepts.values()), extrapolate=extra)
    
    
    def get_t(self, P1, Y1) -> float:
        '''Get temperature from pressure and property value'''
        self.check_range(P1=P1, Y1=Y1)
        func = self.func_t(P1)
        solns = self.combine(func.solve(Y1))
        if len(solns)!=1:
            print(solns)
            raise ValueError("No Unique Solution Found")
        return solns[0]
    
    def combine(self, vals: np.array) -> list:
        '''Remove duplicate and NaN solutions'''
        unique = []
        for i in range(len(vals)):
            for j in range(i,len(vals)):
                if i!=j and (abs(vals[i] - vals[j]) < self.tol):
                    break
                if np.isnan(vals[i]):
                    break
            else:
                unique.append(vals[i])
        return unique

class FluidProperty:
    '''Representation of a fluid property such as pressure or entropy'''
    def __init__(self, tag:str, name:str, alt:tuple, rev:bool=True, file:str=None, units:str=None, axis=False):
        self.name = name
        self.axis = axis
        self.tag = tag
        self.units = units
        if isinstance(alt, str):
            alt = [alt]
        else:
            alt = list(alt)
        self.altnames = alt
        self.file = file
        self.reversible = rev
        self.known = (axis or (self.file is not None))
        self.callnames = [self.tag.upper(), self.name.upper().replace(" ", "_")] + self.altnames
        
    def read_table(self, fname=None):
        if fname is not None:
            self.file = fname
        if self.file is not None:
            self.known = True
            with open(self.file, mode='r', newline='') as f:
                reader = csv.reader(f)
                data = np.array([[s.replace(',', '') for s in r] for r in reader])
            self.table = Table(self, P=data[0,1:].astype(float), T=data[1:,0].astype(float), data=data[1:,1:].astype(float))

    def knownas(self, string):
        return (string.upper().replace(" ", "_") in self.callnames)

    def __str__(self):
        return f"[{self.tag}] {self.name} ({self.units})"
    
    def printval(self, value):
        return str(self) + ": " + str(value)

class FluidState:
    '''Representation of one particular state of a fluid with a given pressure, temperature, etc.'''
    def __init__(self, fluid, props):
        self.properties = dict()
        self.fluid = fluid # Corresponding fluid with property tables
        self.solve_properties(props) # solve for fluid properties

    def __getattr__(self, __name: str):
        '''overrides default FluidState.attribute behavior for unrecognized attribute names. Allows for parsing of '''
        for p in self.fluid.properties:
            if p.knownas(__name) and p.known:
                return self.properties[p]
    
    def solve_properties(self, props: dict):
        '''Solve for all properties of the fluid'''
        pt_count = 0
        keys = list(props.keys())
        if self.fluid.v in keys: # convert specific volume to density
            props[self.fluid.rho] = 1 / props[self.fluid.volume]
            keys.remove(self.fluid.volume)
            keys.append(self.fluid.rho)

        for pr, val in props.items():
            if pr.known:
                if pr is self.fluid.P:
                    pt_count +=1
                elif pr is self.fluid.T:
                    pt_count +=1
                elif not pr.reversible:
                    raise ValueError(f"Property '{pr}' cannot be used to solve for other properties")
            else:
                raise ValueError(f"Property '{pr}' values do not exist for this fluid")

        # Solve for both pressure and temperature if needed
        if pt_count==0:
            self.double_solve(props[keys[0]], keys[0].table, props[keys[1]], keys[1].table)
        
        # Solve for only one of either pressure or temperature if needed
        if pt_count == 1:
            # Check if temperature was provided
            if self.fluid.T in keys:
                # solve for pressure
                self.properties[self.fluid.T] = props[self.fluid.T]
                other_prop = keys[keys is not self.fluid.T]
                self.properties[self.fluid.P] = other_prop.table.get_p( T1=props[self.fluid.T], Y1=props[other_prop])

            # check if pressure was provided 
            elif self.fluid.P in keys:
                # self for temperature
                self.properties[self.fluid.P] = props[self.fluid.P]
                other_prop = keys[keys is not self.fluid.P]
                self.properties[self.fluid.T] = other_prop.table.get_t( P1=props[self.fluid.P], Y1=props[other_prop])
            else:
                # neither property can be found but one was flagged as being provided in previous section
                raise NameError("Solve Error: Could not locate provided pressure/temperature")
        
        # set value for each property
        for pr in self.fluid.properties:
            if pr in keys:
                self.properties[pr] = props[pr]
            elif pr.known and (not pr.axis):
                try:
                    self.properties[pr] = pr.table.interp(P1=self.P, T1=self.T)
                except KeyError as err:
                    # TODO: FIXME
                    
                    print("FAILED:", err, pr.name, "P,T:", self.P, self.T)
                    #print(pr.table.cs_t)
                    self.properties[pr] = 0
                    #self.properties[pr] = np.nan
    
    def double_solve(self, value1, table1, value2, table2) -> "float, float":
        '''solve for both pressure and temperature as a function of other variables
        Should mostly only be called internally by other fluid state methods'''

        t1_p = table1.t_func_p(value1) # function T(P) given property 1
        t2_p = table2.t_func_p(value2) # function T(P) given property 2
        
        pmax = min(max(t1_p.x), max(t2_p.x)) # upper bound of smaller range
        pmin = max(min(t1_p.x), min(t2_p.x)) # lower bound of smaller range
        x_new = []
        c_new = []
        dbp("P RANGE:", pmin, pmax)
        # iterate over pressure values for each
        dbp("X1:", t1_p.x.shape)
        dbp("C1:", t1_p.c.shape)
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
                dbp("P1:", p1, "P2:", p2)
                
            else:
                continue
            break
    
        c=np.array(c_new).transpose(); x=np.array(x_new)
        dbp(c.shape)
        dbp(x.shape)
        difference = PPoly(c=c, x=x, extrapolate=False)
        pressure_values = table1.combine(difference.solve(0))
        if len(pressure_values)!=1:
            print(pressure_values)
            raise ValueError("No Unique Solution Found")
        else:
            temperature = t1_p(pressure_values[0])
            self.properties[self.fluid.P] = pressure_values[0]
            self.properties[self.fluid.T] = temperature
            #return pressure_values[0], temperature
                
    # control printout behavior of fluid state
    def __str__(self):

        pstr = ["Fluid state with properties:"] + [f"{p}: {self.properties[p]}" for p in self.properties.keys()]
        return "\n\t".join(pstr)

class Fluid:

    '''Representation of a fluid medium. Contains data for the fluid across different states of varying pressure, temperature, etc.'''
    def __init__(self, name: str, prop_tables:'{str:str}'):
        self.name = name
        self.properties = [FluidProperty(tag='p', name='Pressure', alt=('PRESS'), rev=False, axis=True),    # Pressure
            FluidProperty(tag='t', name='Temperature', alt=('TEMP'), rev=False, axis=True),                 # Temperature
            FluidProperty(tag='h', name='Enthalpy', alt=('ENTH'), rev=True),                                # Enthalpy
            FluidProperty(tag='u', name='Internal Energy', alt=(), rev=True),                               # Internal Energy
            FluidProperty(tag='v', name='Specific Volume', alt=('VOLUME', 'VOL'), rev=True),                # Specific volume
            FluidProperty(tag='rho', name='Density', alt=('R'), rev=True),                                  # Density
            FluidProperty(tag='s', name='Entropy', alt=('S', 'ENTROPY'), rev=True),                         # Entropy
            FluidProperty(tag='cp', name='Specific Heat', alt=('C'), rev=False),                            # Specific heat
            FluidProperty(tag='gamma', name='Specific Heat Ratio', alt=('Y'), rev=False),                   # gamma, sp. heat ratio
            FluidProperty(tag='nu', name='Kinematic Viscosity', alt=('VISC', 'VISCOSITY'), rev=False),      # viscosity (kinematic)
            FluidProperty(tag='a', name='Speed of Sound', alt=(), rev=False),                               # speed of sound
            FluidProperty(tag='gibbs', name='Gibbs Free Energy', alt=('GFE'), rev=False),                   # gibbs free energy
            FluidProperty(tag='k', name='Thermal Conductivity', alt=("CONDUCTIVITY"), rev=False),                         # thermal conductivity
            FluidProperty(tag='m', name='Molecular Mass', alt=("MOLAR_MASS"), rev=False)]                   # molecular mass
        for prop, fname in prop_tables.items():
            for p in self.properties:
                if p.knownas(prop):
                    #print(p)
                    p.read_table(fname)
        
    
    # Control printout behavior
    def __str__(self):
        return f"Fluid object '{self.name}'"
    
    def __call__(self, **kwargs: float) -> FluidState:
        return self.state(**kwargs)


    def __getattr__(self, __name: str):
        for p in self.properties:
            if p.knownas(__name):
                return p

    
    def state(self, **kwargs: float) -> FluidState:
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
        return FluidState(self, input_props)

class Process:
    def __init__(self, state1, **kwargs):
        self.isentropic = None
        self.adiabatic = None
        self.isenthalpic = None
        self.q = 0
        self.w = 0
        
        # TODO: Define process which goes from one state to another
        pass


#dir = "CEA Property Tables"
dir = "Updated"
prop_files = dict()
for file in os.listdir(dir):
    tag = file.split(".")[0]
    prop_files[tag] = os.path.join(dir, file)

