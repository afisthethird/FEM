# Libraries
import sympy as sp
from dataclasses import dataclass
from enum import IntEnum
# Scripts
from Code.types import *
from Code.symbolic.space import Space
from Code.symbolic.math import SymbolicMathematicalObject
from Code.utilities.auxilary import make_callable, extract_base_types


'''
Script-specific typing setup
'''

LessThanBoundaryRelationalType : TypeAlias = Union[
    sp.Le, 
    sp.Lt
    ]
GreaterThanBoundaryRelationalType : TypeAlias = Union[
    sp.Ge, 
    sp.Gt
    ]
BoundaryRelationalType : TypeAlias = Union[
    LessThanBoundaryRelationalType, 
    GreaterThanBoundaryRelationalType
    ]
RUNTIME_LESS_THAN_BOUNDARY_RELATIONAL_TYPE = extract_base_types(LessThanBoundaryRelationalType)
RUNTIME_GREATER_THAN_BOUNDARY_RELATIONAL_TYPE = extract_base_types(GreaterThanBoundaryRelationalType)
RUNTIME_BOUNDARY_RELATIONAL_TYPE = extract_base_types(BoundaryRelationalType)

CoordinateType : TypeAlias = NumericRowVectorValueType

# Bookkeeping "class" to make coordinate location returns more clear
class CoordinateLocated(IntEnum):
    OUTSIDE = 0
    INSIDE = 1
    ON = 2


'''
Mathematical Region Objects
'''

class BoundarySurface:
    """
    Essentially a wrapper for a SymPy `Expr` dictating the boundary in space, as a function of that `Space`'s dimensional symbols.
    Takes a SymPy relational object as input, and converts it to the standardized form for sign checking.

    Ex: `sp.Expr` `g`; for (vector) coordinate `v` such that the callable (lambdified) `g(v)` ≤ 0 ⟹ `v` is inside that boundary.
    """

    def __init__(
        self,
        name : NameType,
        eq : BoundaryRelationalType,
        host_spce : Space,
        cndtn : SymbolicMathematicalObject,
        tol : NumericDecimalValueType = 1e-10
        ):

        self.name = name

        self.eq = eq # Kept for display purposes
        self.host_spce = host_spce
        self.fn = make_callable(
            self.host_spce.dims_syms(), 
            self.convert_equation_to_expression(eq)
            )
        self.cndtn = cndtn
        self.tol = tol

    # Converts the user's input equation describing the boundary to a standard-form expression:
    # LHS ~ RHS ⟼ f(LHS, RHS) ≤ 0
    @staticmethod
    def convert_equation_to_expression(
        eq : BoundaryRelationalType
        ) -> sp.Expr:
        
        eq_lhs = eq.lhs
        eq_rhs = eq.rhs
        expr = None
        if isinstance(eq, RUNTIME_LESS_THAN_BOUNDARY_RELATIONAL_TYPE):
            expr = eq_lhs - eq_rhs
        elif isinstance(eq, RUNTIME_GREATER_THAN_BOUNDARY_RELATIONAL_TYPE):
            expr = eq_rhs - eq_lhs
        else:
            raise ValueError(f"Provided boundary-defining equation {eq} " +
                f"is not one of the following supported relational types: {RUNTIME_BOUNDARY_RELATIONAL_TYPE}"
                )

        return expr

    # 0 for outside, 1 for inside, 2 for on
    def contains(
        self,
        crd : CoordinateType
        ) -> NumericIntegerValueType:

        if len(crd) != len(self.host_spce):
            raise ValueError(
                f"Dimensionality of coordinate {crd} " +
                f"does not match dimensionality of boundary \"{self.name}\" " +
                f"with host space \"{self.host_spce.name}\"."
                )

        fn_result = self.fn(*crd)
        if fn_result >= self.tol:
            result = CoordinateLocated.OUTSIDE
        elif fn_result <= -self.tol:
            result = CoordinateLocated.INSIDE
        else:
            result = CoordinateLocated.ON
        # Typecast to int to make parsing & printing cleaner
        result = int(result)

        return result


# Closed domain defined by boundaries
@dataclass
class Domain:

    name : NameType
    bdrys : List[BoundarySurface]
    host_spce : Space

    def contains(
        self, 
        crd : CoordinateType
        ) -> List[NumericIntegerValueType]:

        if len(crd) != len(self.host_spce):
            raise ValueError(
                f"Dimensionality of coordinate {crd} " +
                f"does not match dimensionality of domain \"{self.name}\" " +
                f"with host space \"{self.host_spce.name}\"."
                )
        
        result = [
            curr_bdry.contains(crd)
            for curr_bdry in self.bdrys
            ]
        
        return result


'''
Defaults
'''

