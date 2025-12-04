# Libraries
from abc import ABC, abstractmethod
# Scripts
from ...types import *


'''
Script-specific typing setup
'''

if TYPE_CHECKING:
    from ..base import CompactTopologicalSpace


'''
ConnectivityStructure superclasses
'''

class ConnectivityStructure(ABC):
    DIMALTY : NumericIntegerValueType
    NUM_OF_CMPNTS : NumericIntegerValueType
    CMPNTS_CNCTVTY : Tuple[Tuple[IndexType]] = None
    CMPNT_TYPES : Tuple[Type["ConnectivityStructure"]] = None

    def __init__(
        self,
        host_spce : "CompactTopologicalSpace",
        nums_of_cmpnts_gend_per_host_spce_dim : Dict[IndexType, NumericIntegerValueType] = None,
        idx : IndexType = None,
        cmpnt_subs : Dict[IndexType, "ConnectivityStructure"] = None
        ):

        # Host space
        self.host_spce = host_spce

        # Indexing handling and defaults
        # If conflicting arguments are provided
        if (idx is not None) and (nums_of_cmpnts_gend_per_host_spce_dim is not None):
            raise ValueError("Overlapping indexing values provided.")
        # If nums_of_cmpnts_gend_per_topo_dimalty was provided
        elif nums_of_cmpnts_gend_per_host_spce_dim is not None:
            self.idx = nums_of_cmpnts_gend_per_host_spce_dim[self.DIMALTY]
            nums_of_cmpnts_gend_per_host_spce_dim[self.DIMALTY] += 1
        # If idx was provided (should only occur for modification of existing Topologies)
        elif idx is not None:
            self.idx = idx
        # If only host_spce was provided; no indexing provided yet, use default
        else:
            nums_of_cmpnts_gend_per_host_spce_dim = {
                curr_dim_idx: 0 for curr_dim_idx in range(self.host_spce.DIMALTY+1)
                }
            self.idx = nums_of_cmpnts_gend_per_host_spce_dim[self.DIMALTY]
            nums_of_cmpnts_gend_per_host_spce_dim[self.DIMALTY] += 1

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
                        host_spce = self.host_spce,
                        nums_of_cmpnts_gend_per_host_spce_dim = nums_of_cmpnts_gend_per_host_spce_dim,
                        cmpnt_subs = subs_for_curr_cmpnt
                        )
                # Component definition, for cases of no subcomponents
                else:
                    curr_cmpnt = curr_cmpnt_type(
                        host_spce = self.host_spce,
                        nums_of_cmpnts_gend_per_host_spce_dim = nums_of_cmpnts_gend_per_host_spce_dim,
                        )
            self.cmpnts.append(curr_cmpnt)

    def __eq__(self, other) -> bool:
        if not isinstance(other, type(self)):
            return False
        return all(other_curr_cmpnt in self.cmpnts for other_curr_cmpnt in other.cmpnts)

    def __getitem__(self, key: IndexType) -> "ConnectivityStructure":
        return self.cmpnts[key]

    @abstractmethod
    def __neg__(self) -> "ConnectivityStructure":
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

# 0D : Points
class PointStructure(ConnectivityStructure, ABC):
    DIMALTY = 0
    NUM_OF_CMPNTS = 0

    # There is no sense of inversion for a 0D topology
    def __neg__(self) -> "PointStructure":
        return self

# 1D : Curves
class CurveStructure(ConnectivityStructure, ABC):
    DIMALTY = 0
    CMPNT_TYPES : Tuple[Type[PointStructure]]

    # Negating/inverting a curve is equivalent to reversing its direction
    def __neg__(self) -> "CurveStructure":
        cls = type(self)
        return cls(
            host_spce = self.host_spce,
            idx = self.idx,
            cmpnt_subs = dict(enumerate(self.cmpnts[::-1]))
            )

# 2D : Surfaces
class SurfaceStructure(ConnectivityStructure, ABC):
    DIMALTY = 2
    CMPNT_TYPES : Tuple[Type[CurveStructure]]

    # There is currently no implementation of 4D Topologies, so Surface does not need a negation/inversion operator defined
    def __neg__(self):
        raise NotImplementedError

# 3D : Volumes
class VolumeStructure(ConnectivityStructure, ABC):
    DIMALTY = 3
    CMPNT_TYPES : Tuple[Type[SurfaceStructure]]

    # There is currently no implementation of 4D Topologies, so SurfaceTopology does not need a negation/inversion operator defined
    def __neg__(self):
        raise NotImplementedError

