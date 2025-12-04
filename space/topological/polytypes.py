# Libraries
from abc import ABC
from dataclasses import dataclass
# Scripts
from ...types import *
from .connectivity import ConnectivityStructure, PointStructure, CurveStructure, SurfaceStructure, VolumeStructure


'''
Superclass
'''

class PolytypeStructure(ConnectivityStructure, ABC):
    CMPNT_TYPES : Tuple[Type["PolytypeStructure"]] = None


'''
0D : Atomic
'''

# 0D polytype structure class (Vertices are atomic, and thus not unique in structure)
@dataclass(init=False, repr=False)
class Vertex(PolytypeStructure, PointStructure):
    pass


'''
1D : Composed of 0D
'''

# 1D polytype structure class (Edges are atomic, and thus not unique in structure)
@dataclass(init=False, repr=False)
class Edge(PolytypeStructure, CurveStructure):
    NUM_OF_CMPNTS = 2
    CMPNT_TYPES = (
        Vertex,
        Vertex
        )


'''
2D : Composed of 1D
'''

# 2D polytype structure superclass
class Face(PolytypeStructure, SurfaceStructure, ABC):
    CMPNT_TYPES : Tuple[Type[Edge]]

@dataclass(init=False, repr=False)
class Quadrilateral(Face):
    NUM_OF_CMPNTS = 4
    CMPNT_TYPES = (
        Edge,
        Edge,
        Edge,
        Edge
        )
    CMPNTS_CNCTVTY = (
        (3, 1),
        (0, 2),
        (1, 3),
        (2, 0)
        )

@dataclass(init=False, repr=False)
class Triangle(Face):
    NUM_OF_CMPNTS = 3
    CMPNT_TYPES = (
        Edge,
        Edge,
        Edge
        )
    CMPNTS_CNCTVTY = (
        (2, 1),
        (0, 2),
        (1, 0)
        )


'''
3D : Composed of 2D
'''

# 3D polytype structure superclass
class Polyhedron(PolytypeStructure, VolumeStructure, ABC):
    CMPNT_TYPES : Tuple[Type[Face]]

@dataclass(init=False, repr=False)
class TriangularPrism(Polyhedron):
    NUM_OF_CMPNTS = 5
    CMPNTS_CNCTVTY = (
        (1, 2, 3, 4),
        (0, 4, 2),
        (0, 1, 4, 3),
        (0, 2, 4),
        (0, 3, 2, 1)
        )
    CMPNT_TYPES = (
        Quadrilateral,
        Triangle,
        Quadrilateral,
        Triangle,
        Quadrilateral
        )

@dataclass(init=False, repr=False)
class Pyramid(Polyhedron):
    NUM_OF_CMPNTS = 5
    CMPNTS_CNCTVTY = (
        (1, 2, 3, 4),
        (0, 4, 2),
        (0, 1, 3),
        (0, 2, 4),
        (0, 3, 1)
        )
    CMPNT_TYPES = (
        Quadrilateral,
        Triangle,
        Triangle,
        Triangle,
        Triangle
        )

