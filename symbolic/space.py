# Libraries
import sympy as sp
from abc import ABC
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


# Currently a dataclass but may be extended to non-Cartesian coordinate systems in the future
@dataclass
class Space(ABC):
    name : NameType
    dims : Tuple[Dimension, ...]
    
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

@dataclass(repr=False)
class PhysicalSpace(Space):
    ref_spce : "ReferenceSpace" = None

    def create_reference_space(
        self,
        name : NameType = None,
        dims : Tuple[Dimension, ...] = None    
        ):

        # Don't overwrite
        if self.ref_spce is not None:
            return
        # Default reference space naming scheme
        if name is None:
            name = f"{self.name}_ref"
        if dims is None:
            dims = tuple([
                Dimension(sp.Symbol(f"{curr_dim}_ref"))
                for curr_dim in self.dims
                ])
        # Creation and storage
        self.ref_spce = ReferenceSpace(
            name = name,
            dims = dims,
            phys_spce = self
            )

@dataclass(repr=False)
class ReferenceSpace(Space):
    phys_spce : PhysicalSpace = None


'''
Default Cartesian coordinate spaces
'''

x_dim = Dimension(sp.Symbol('x'))
x_ref_dim = Dimension(sp.Symbol('ξ'))

R1 = PhysicalSpace(
    name = 'R1', 
    dims = (x_dim,)
    )
R1.create_reference_space(
    name = 'R1_ref',
    dims = (x_ref_dim,)
    )


y_dim = Dimension(sp.Symbol('y'))
y_ref_dim = Dimension(sp.Symbol('η'))

R2 = PhysicalSpace(
    name = 'R2',
    dims = (x_dim, y_dim)
    )
R2.create_reference_space(
    name = 'R2_ref', 
    dims = (x_ref_dim, y_ref_dim)
    )


z_dim = Dimension(sp.Symbol('z'))
z_ref_dim = Dimension(sp.Symbol('ζ'))

R3 = PhysicalSpace(
    name = 'R3', 
    dims = (x_dim, y_dim, z_dim)
    )
R3.create_reference_space(
    name = 'R3_ref', 
    dims = (x_ref_dim, y_ref_dim, z_ref_dim)
    )

