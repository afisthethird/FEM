# Libraries
from abc import ABC, abstractmethod
from dataclasses import dataclass
# Scripts
from Code.types import *


'''
Topology superclass
'''

class Topology(ABC):
    NUM_OF_DIMS : NumericIntegerValueType
    NUM_OF_CMPNTS : NumericIntegerValueType
    CMPNT_TYPES : Tuple["Topology"] = None
    CMPNTS_CNCTVTY : Tuple[Tuple[IndexType]] = None

    def __init__(
        self,
        idx : IndexType = None,
        nums_of_topos_gend_per_dim : Dict[IndexType, NumericIntegerValueType] = None,
        cmpnt_subs : Dict[IndexType, "Topology"] = None
        ):

        # Indexing handling and defaults
        if (idx is not None) and (nums_of_topos_gend_per_dim is not None):
            raise ValueError("Overlapping indexing values provided.")
        elif idx is not None:
            self.idx = idx
        else:
            # No indexing provided yet, use default
            if nums_of_topos_gend_per_dim is None:
                nums_of_topos_gend_per_dim = {
                    curr_dim_idx: 0 for curr_dim_idx in range(self.NUM_OF_DIMS+1)
                    }
            self.idx = nums_of_topos_gend_per_dim[self.NUM_OF_DIMS]
            nums_of_topos_gend_per_dim[self.NUM_OF_DIMS] += 1

        # Component generation
        self.cmpnts = []
        for curr_cmpnt_idx in range(self.NUM_OF_CMPNTS):
            curr_cmpnt = None
            # Check provided substitutions
            if (cmpnt_subs is not None) and (curr_cmpnt_idx in cmpnt_subs):
                curr_cmpnt = cmpnt_subs[curr_cmpnt_idx]
            else:
                curr_cmpnt_type = self.CMPNT_TYPES[curr_cmpnt_idx]
                # Subcomponent substitution, if applicable
                if curr_cmpnt_type.NUM_OF_CMPNTS > 0:
                    subs_for_curr_cmpnt = {}
                    # Use connectivity vectors to identify shared subcomponents
                    for curr_cnctd_cmpnt_idx in self.CMPNTS_CNCTVTY[curr_cmpnt_idx]:
                        if curr_cnctd_cmpnt_idx < curr_cmpnt_idx:
                            subs_for_curr_cmpnt[
                                self.CMPNTS_CNCTVTY[curr_cmpnt_idx].index(curr_cnctd_cmpnt_idx)
                                ] = -self.cmpnts[curr_cnctd_cmpnt_idx][
                                    self.CMPNTS_CNCTVTY[curr_cnctd_cmpnt_idx].index(curr_cmpnt_idx)
                                    ]
                    # Current component definition when subcomponents exist
                    curr_cmpnt = curr_cmpnt_type(
                        nums_of_topos_gend_per_dim = nums_of_topos_gend_per_dim,
                        cmpnt_subs = subs_for_curr_cmpnt
                        )
                # Component definition, for cases of no subcomponents
                else:
                    curr_cmpnt = curr_cmpnt_type(
                        nums_of_topos_gend_per_dim = nums_of_topos_gend_per_dim,
                        )
            self.cmpnts.append(curr_cmpnt)

    def __contains__(self, item) -> bool:
        if isinstance(item, Topology):
            if item.NUM_OF_DIMS == (self.NUM_OF_DIMS - 1):
                return item in self.cmpnts
            elif item.NUM_OF_DIMS < (self.NUM_OF_DIMS - 1):
                return any(item in curr_cmpnt for curr_cmpnt in self.cmpnts)
        return False

    def __eq__(self, other) -> bool:
        if not isinstance(other, type(self)):
            return False
        return all(other_curr_cmpnt in self.cmpnts for other_curr_cmpnt in other.cmpnts)

    def __getitem__(self, key: IndexType) -> "Topology":
        return self.cmpnts[key]

    @abstractmethod
    def __neg__(self) -> "Topology":
        pass
    
    def __repr__(self):
        repr = f"{type(self).__name__}{self.idx}"
        if self.cmpnts:
            repr += f"({", ".join([
                curr_cmpnt.__repr__() 
                for curr_cmpnt 
                in self.cmpnts
                ])})"
        return repr
        

'''
0D Topologies : Points
'''

class Point(Topology, ABC):
    NUM_OF_DIMS = 0
    NUM_OF_CMPNTS = 0
    
    # There is no sense of inversion for a 0D topology
    def __neg__(self) -> "Point":
        return self

@dataclass
class Vertex(Point):
    pass


'''
1D Topologies : Curves
'''

class Curve(Topology, ABC):
    NUM_OF_DIMS = 1
    CMPNT_TYPES : Tuple[Type[Point]]
    
    # Negating/inverting a curve is equivalent to reversing its direction
    def __neg__(self) -> "Curve":
        cls = type(self)
        return cls(
            idx = self.idx,
            cmpnt_subs = dict(enumerate(self.cmpnts[::-1]))
            )

# Planar Curves
@dataclass
class Edge(Curve):
    NUM_OF_CMPNTS = 2
    CMPNT_TYPES = (
        Vertex,
        Vertex
        )

    
'''
2D Topologies : Surfaces
'''

class Surface(Topology, ABC):
    NUM_OF_DIMS = 2
    CMPNT_TYPES : Tuple[Type[Curve]]

    # There is currently no implementation of 4D Topologies, so Surface does not need a negation/inversion operator defined
    def __neg__(self):
        raise NotImplementedError

# Planar Surfaces
@dataclass
class Face(Surface, ABC):
    CMPNT_TYPES : Tuple[Type[Edge]]

@dataclass
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

@dataclass
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
3D Topologies : Volumes
'''

class Volume(Topology, ABC):
    NUM_OF_DIMS = 3
    CMPNT_TYPES : Tuple[Type[Surface]]

    # There is currently no implementation of 4D Topologies, so Surface does not need a negation/inversion operator defined
    def __neg__(self):
        raise NotImplementedError

@dataclass
class TriangularPrism(Volume):
    NUM_OF_CMPNTS = 5
    CMPNT_TYPES = (
        Quadrilateral,
        Triangle,
        Quadrilateral,
        Triangle,
        Quadrilateral
        )
    CMPNTS_CNCTVTY = (
        (1, 2, 3, 4),
        (0, 4, 2),
        (0, 1, 4, 3),
        (0, 2, 4),
        (0, 3, 2, 1)
        )

@dataclass
class Pyramid(Volume):
    NUM_OF_CMPNTS = 5
    CMPNT_TYPES = (
        Quadrilateral,
        Triangle,
        Triangle,
        Triangle,
        Triangle
        )
    CMPNTS_CNCTVTY = (
        (1, 2, 3, 4),
        (0, 4, 2),
        (0, 1, 3),
        (0, 2, 4),
        (0, 3, 1)
        )

