# Libraries
import sympy as sp
from abc import ABC
from dataclasses import dataclass
# Scripts
from Code.types import *


# Individual "coordinates" are mathematically nothing more than consistently-labeled placeholders- that is, a "coordinate" identifies a particular “slot” corresponding to the basis vector scaled by its supplied coordinate vector's concrete argument.
@dataclass(frozen=True)
class Coordinate:
    sym : sp.Symbol
    
    def __hash__(self):
        return hash(self.sym)

    def __eq__(self, other):
        if isinstance(other, Coordinate):
            if other.sym == self.sym:
                return True
        return False

    def __repr__(self):
        return f"{self.sym}"

class CoordinateVector(Vector):
    host_spce : "CoordinateSpace"
    cmpnts : Dict[Coordinate, SymbolicScalarValueType]

    # Input checking
    def __post_init__(self):
        host_spce_crds = self.host_spce.crds
        cmpnts_crds = set(self.cmpnts.keys())
        if host_spce_crds != cmpnts_crds:
            raise ValueError(f"Coordinates provided ({cmpnts_crds}) do not match host space coordinates ({host_spce_crds}).")

    def __getitem__(self, key: Coordinate) -> SymbolicScalarValueType:
        return self.cmpnts[key]
    def __setitem__(self, key: Coordinate, val: SymbolicScalarValueType):
        self.cmpnts[key] = val

    @property
    def crds(self) -> Set[Coordinate]:
        return self.host_spce.crds
    @property
    def crds_syms(self) -> Dict[Coordinate, sp.Symbol]:
        return {curr_crd: curr_crd.sym for curr_crd in self.crds}

    def __add__(self, other: "CoordinateVector") -> "CoordinateVector":
        # Input checking
        if other not in self.host_spce:
            raise ValueError(f"Could not sum {self} of space {self.host_spce} with {other} of space {other.host_spce}.")
        return CoordinateVector(
            host_spce = self.host_spce,
            cmpnts = {
                curr_crd: self.cmpnts[curr_crd] + other.cmpnts[curr_crd]
                for curr_crd in self.crds
                }
            )
    
    def __mul__(self, other: SymbolicScalarValueType) -> "CoordinateVector":
        return CoordinateVector(
            host_spce = self.host_spce,
            cmpnts = {
                curr_crd: other * self.cmpnts[curr_crd]
                for curr_crd in self.crds
                }
            )

    def __repr__(self):
        return ", ".join(
            f"{curr_crd}:{self[curr_crd]}" 
            for curr_crd in self.crds
            )

@dataclass(frozen=True)
class CoordinateSpace(VectorSpace, ABC):
    crds : Set[Coordinate]

    # Input checking
    def __post_init__(self):
        crds_dimalty = len(self.crds)
        if self.dimalty != crds_dimalty:
            raise ValueError(f"Number of coordinates provided ({crds_dimalty}) does not match dimensionality of host space ({self.dimalty}).")

    def vector(
        self,
        cmpnts : Dict[Coordinate, SymbolicScalarValueType]
        ):

        

    # Cartesian coordinate space is used as the identity CoordinateSpace
    @abstractmethod
    def identity_isomorphism(
        self, 
        crd_vec: "CoordinateVector"
        ) -> "CoordinateVector":

        pass

    @cached_property
    def basis(self) -> Dict[Coordinate, "CoordinateVector"]:


        





# class 𝑋1(CartesianCoordinateSpace):






quit()




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
Default spaces: 
'''


x_crd = Coordinate(sp.Symbol('x'))
y_crd = Coordinate(sp.Symbol('y'))
z_crd = Coordinate(sp.Symbol('z'))
r_crd = Coordinate(sp.Symbol('r'))
θ_crd = Coordinate(sp.Symbol('θ'))
ɸ_crd = Coordinate(sp.Symbol('ɸ'))
ϕ_crd = Coordinate(sp.Symbol('ϕ'))
ξ_crd = Coordinate(sp.Symbol('ξ'))
η_crd = Coordinate(sp.Symbol('η'))
ζ_crd = Coordinate(sp.Symbol('ζ'))





R1_dim = Dimension(sp.Symbol('x'))
R1_ref_dim = Dimension(sp.Symbol('ξ'))

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

