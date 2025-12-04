# Libraries
from abc import ABC, abstractmethod
from dataclasses import dataclass
# Scripts
from Code.types import *
from Code.space.space import Space, Space0D, Space1D, Space2D, Space3D
from Code.space.topological import TopologicalSpace, NoncompactTopologicalSpace, CompactTopologicalSpace


# Goal: order-based node insertion into Topologies
# Goal: abstract/templated assembly of CoordinateSpaces into ManifoldSpaces


@dataclass(frozen=True)
class ManifoldSpace(Space, ABC):
    smooth : bool
    pass