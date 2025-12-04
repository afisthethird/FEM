# Libraries
from abc import ABC
# Scripts
from Code.types import *
from Code.space.base import Space, Space0D, Space1D, Space2D, Space3D
from Code.space.topological.connectivity import ConnectivityStructure, PointStructure, CurveStructure, SurfaceStructure, VolumeStructure


'''
TopologicalSpace superclasses
'''

class TopologicalSpace(Space, ABC):
    CMPCT : bool

class NoncompactTopologicalSpace(TopologicalSpace, ABC):
    CMPCT = False

class CompactTopologicalSpace(TopologicalSpace, ABC):
    CMPCT = True

    def __init__(self, cnctvty_struct_type: Type["ConnectivityStructure"]):
        # Connectivity information implictly encodes topology set definition
        self.cnctvty_struct = cnctvty_struct_type(host_spce = self)


'''
TopologicalSpaces by dimensionality
'''

# 0D host TopologicalSpace
class Point(Space0D, CompactTopologicalSpace):
    def __init__(self, cnctvty_struct_type: Type["PointStructure"]):
        super().__init__(cnctvty_struct_type = cnctvty_struct_type)

# 1D host TopologicalSpace
class Curve(Space1D, CompactTopologicalSpace):
    def __init__(self, cnctvty_struct_type: Type["CurveStructure"]):
        super().__init__(cnctvty_struct_type = cnctvty_struct_type)

# 2D host TopologicalSpace
class Surface(Space2D, CompactTopologicalSpace):
    def __init__(self, cnctvty_struct_type: Type["SurfaceStructure"]):
        super().__init__(cnctvty_struct_type = cnctvty_struct_type)

# 3D host TopologicalSpace
class Volume(Space3D, CompactTopologicalSpace):
    def __init__(self, cnctvty_struct_type: Type["VolumeStructure"]):
        super().__init__(cnctvty_struct_type = cnctvty_struct_type)

