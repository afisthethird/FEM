
@dataclass
class Vector(ABC):
    host_spce : "VectorSpace"
    cmpnts : Dict[TypeVar, SymbolicScalarValueType]

@dataclass(frozen=True)
class VectorSpace(ManifoldSpace, ABC):
    pass

@dataclass(frozen=True)
class R0(VectorSpace, Point):
    pass
@dataclass(frozen=True)
class R1(VectorSpace, Curve):
    pass
@dataclass(frozen=True)
class R2(VectorSpace, Surface):
    pass
@dataclass(frozen=True)
class R3(VectorSpace, Volume):
    pass

