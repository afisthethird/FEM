# Libraries
import sympy as sp
import numpy as np
from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import Enum
# Scripts
from Code.types import *
from Code.symbolic.space import Space, Dimension


'''
Script-specific typing setup
'''

@dataclass
class Coordinate:

    name : NameType
    dims : Tuple[Dimension, ...]
    vals : Tuple[ScalarValueType]

    def __getitem__(
        self, 
        key : Union[IndexType, Dimension]
        ) -> NumericValueType:

        result_idx = None
        if isinstance(key, IndexType):
            result = self.vals[key]
        elif isinstance(key, Dimension):
            result_idx = self.dims.index(key)
        
        result = None
        if result_idx is not None:
            result = self.vals[result_idx]
        
        return result
    
    def __len__(self):

        return len(self.dims)

@dataclass
class Node:

    name : NameType
    phys_crd : Coordinate
    ref_crd : Coordinate



# Utility parent class for standardized definition of element topologies
@dataclass(frozen=True)
class ElementTopology(ABC):

    num_of_verts : NumericIntegerValueType
    edges_verts_cnctvty : np.ndarray[NumericIntegerValueType, Literal["(total edges, 2)"]]
    faces_verts_cnctvty : np.ndarray[NumericIntegerValueType, Literal["(total faces, vertices per face)"]]
    
    # Node ordering convention: [VERTICES, EDGES, FACES, VOLUME]
    def __init__(
        self,
        order : NumericIntegerValueType,
        edges_orders : Annotated[Tuple[NumericIntegerValueType, ...], Literal["(total edges)"]] = None,
        faces_orders : Annotated[Tuple[NumericIntegerValueType, ...], Literal["(total faces)"]] = None,
        vol_order : NumericIntegerValueType = None
        ):

        self.num_of_nds_per_edge = self.generate_number_of_nodes_per_edge()
        self.vert_nds_idxs = self.generate_vertex_nodes_indices()
        self.edge_nds_cnctvty_idxs = self.generate_edge_nodes_connectivity_indices()
        self.face_nds_cnctvty_idxs = self.generate_face_nodes_connectivity_indices()

    @abstractmethod
    def compute_number_of_nodes(self):
        pass
        
    @abstractmethod
    def generate_vertex_nodes_indices(self):
        pass

    @abstractmethod
    def generate_edge_nodes_connectivity_indices(self):
        pass

    @abstractmethod
    def generate_face_nodes_connectivity_indices(self):
        pass


class ElementSelection(Enum):
    LINE = "Line"


def create_element(
    type : ElementSelection,
    order : int    
    ):

    cls = hi