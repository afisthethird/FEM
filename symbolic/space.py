# Libraries
import sympy as sp
from dataclasses import dataclass
# Scripts
from Code.types import *


@dataclass
class Dimension:
    
    sym : sp.Symbol

    def __repr__(self):

        return f"{self.sym}"

    def __hash__(self):

        return hash(self.sym)

    def __eq__(self, other):

        result = False
        if isinstance(other, Dimension):
            if other.sym == self.sym:
                result = True

        return result


class Space:

    def __init__(
        self,
        name : NameType,
        dims : Tuple[Dimension, ...]
        ):
        
        self.name = name
        self.dims = dims
    
    def __repr__(self):

        return f"{self.name}: {self.dims}"
    
    def __contains__(
        self,
        dim : Dimension
        ) -> bool:

        result = False
        if dim in self.dims:
            result = True
        
        return result

    def __getitem__(
        self, 
        key : IndexType
        ) -> Dimension:

        return self.dims[key]
    
    def __len__(self):

        return len(self.dims)
    
    def dims_syms(self) -> Tuple[sp.Symbol, ...]:

        return tuple([curr_dim.sym for curr_dim in self.dims])


class ReferenceSpace(Space):

    def __init__(
        self,
        name : NameType,
        dims : Tuple[Dimension, ...],
        host_spce : Space
        ):

        super().__init__(name, dims)
        self.host_spce = host_spce


'''
Default Cartesian coordinate spaces
'''

x_dim = Dimension(sp.Symbol('x'))
x_ref_dim = Dimension(sp.Symbol('ξ'))

R1 = Space(
    name = 'R1', 
    dims = (x_dim,)
    )
R1_ref = ReferenceSpace(
    name = 'R1_ref', 
    dims = (x_ref_dim,),
    host_spce = R1
    )


y_dim = Dimension(sp.Symbol('y'))
y_ref_dim = Dimension(sp.Symbol('η'))

R2 = Space(
    name = 'R2',
    dims = (x_dim, y_dim)
    )
R2_ref = ReferenceSpace(
    name = 'R2_ref', 
    dims = (x_ref_dim, y_ref_dim),
    host_spce = R2
    )


z_dim = Dimension(sp.Symbol('z'))
z_ref_dim = Dimension(sp.Symbol('ζ'))

R3 = Space(
    name = 'R3', 
    dims = (x_dim, y_dim, z_dim)
    )
R3_ref = ReferenceSpace(
    name = 'R3_ref', 
    dims = (x_ref_dim, y_ref_dim, z_ref_dim),
    host_spce = R3
    )

