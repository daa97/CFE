class Material:
    def __init__(self, **kwargs):
        Material.properties = ("name", "density", "E", "porosity")
        self.set(**kwargs)


    def set(self, **kwargs):
        for key in self.properties:
            if key in kwargs.keys():
                self.__setattr__(key,kwargs[key])
            else:
                self.__setattr__(key, None)


class Part(Material):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        Part.properties += ("material", "volume", "mass")
        self.set(**kwargs)


class Cylinder(Part):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        Cylinder.properties += ("IR", "ID", "OR", "OD", "length")
        self.set(**kwargs)
        print(self.__dict__)

PM = Cylinder(material="SiC", IR=.045, OR=.049, length=0.84)
CASE = Cylinder(material="304 SS", IR=.051, OR=.056, length=0.84)

