# Libraries
import sympy as sp
from abc import ABC, abstractmethod
from dataclasses import dataclass
# Scripts
from Code.types import *
from Code.symbolic.space import Space


'''
Script-specific typing setup
'''
# Supported automatic composition types
UncastExpressionType : TypeAlias = SymbolicValueType
UncastArgumentType   : TypeAlias = NameType
SymbolicMathematicalObjectInputType : TypeAlias = Union[
    UncastExpressionType, 
    UncastArgumentType, 
    "SymbolicMathematicalObject"
    ]
OperandInputType = SymbolicMathematicalObjectInputType

# Concrete runtime typing for raw operand input casting
RUNTIME_UNCAST_EXPRESSION_TYPE = extract_base_types(UncastExpressionType)
RUNTIME_UNCAST_ARGUMENT_TYPE = extract_base_types(UncastArgumentType)


'''
Abstract Mathematical Object Superclasses
'''

# Abstract superclass for automatic composition overhead
class SymbolicMathematicalObject(ABC):

    # Automatic composition
    def __neg__(self):
        return Negation(self)
    def __add__(self, other: OperandInputType):
        return Sum(self, other)
    def __radd__(self, other: OperandInputType):
        return Sum(other, self)
    def __sub__(self, other: OperandInputType):
        return Difference(self, other)
    def __rsub__(self, other: OperandInputType):
        return Difference(other, self)
    def __mul__(self, other: OperandInputType):
        return ScalarMultiplication(self, other)
    def __rmul__(self, other: OperandInputType):
        return ScalarMultiplication(other, self)
    def __truediv__(self, other: OperandInputType):
        return ScalarDivision(self, other)
    def __rtruediv__(self, other: OperandInputType):
        return ScalarDivision(other, self)
    def dot(self, other: OperandInputType):
        return DotProduct(self, other)
    def __matmul__(self, other: OperandInputType):
        return MatrixProduct(self, other)
    def __rmatmul__(self, other: OperandInputType):
        return MatrixProduct(other, self)
    
    # Method for casting to SymbolicMathematicalObject classes from raw input
    @staticmethod
    def cast(
        inps : List[SymbolicMathematicalObjectInputType]
        ) -> List[Self]:

        cast_inps = []
        for curr_inp in inps:
            curr_cast_inp = None
            if isinstance(curr_inp, RUNTIME_UNCAST_EXPRESSION_TYPE):
                curr_cast_inp = Expression(curr_inp)
            elif isinstance(curr_inp, RUNTIME_UNCAST_ARGUMENT_TYPE):
                curr_cast_inp = Argument(curr_inp)
            elif isinstance(curr_inp, (Expression, Argument, Operator)):
                curr_cast_inp = curr_inp
            else:
                raise TypeError
            cast_inps.append(curr_cast_inp)
            
        return cast_inps
        

'''
Non-Operators
'''

# Wrapper class for representation of expressions embedded in unevaluated operators
@dataclass
class Expression(SymbolicMathematicalObject):

    val : SymbolicValueType

# Wrapper class for represenation of placeholder expressions in unevaluated operators
# Must be supplied for complete evaluation
@dataclass
class Argument(SymbolicMathematicalObject):

    name : NameType

    def __hash__(self):

        return hash(self.name)

    def __eq__(self, other):

        result = False
        if isinstance(other, Argument):
            if other.name == self.name:
                result = True
        elif isinstance(other, RUNTIME_NAME_TYPE):
            if other == self.name:
                result = True

        return result

'''
Operator Superclasses
'''

class Operator(SymbolicMathematicalObject, ABC):
    """
    Base class for unevaluated representation of symbolic mathematical operators.

    This class serves as an abstract template for all concrete operators.
    """

    def __init__(
        self,
        oprnds : List[OperandInputType]
        ):

        # Set cast operands
        self.oprnds = self.cast(oprnds)
    
    def __call__(
        self,
        args : Dict[Union[UncastArgumentType, Argument], SymbolicValueType] = None
        ) -> Tuple[SymbolicValueType]:
        """
        Calling an `Operator` with a list of substitutions for its `Argument` placeholders returns its evaluation.
        The recursive evaluation control logic common to all concrete mathematical operators is contained in this function.
        Specific evaluation logic is reserved for subclass implementations of `evaluate()`.
        
        Note that it is not meaningful to "call" or "evaluate" expressions or arguments, so the `Expression` and `Argument` classes have no such methods.
        Mathematical expressions are by definition already assembled/evaluated, and arguments are nothing but placeholders, so they can not be "called" or "evaluated" in the same sense as mathematical operators.
        """
        
        evald_oprnds = []
        for curr_oprnd in self.oprnds:
            match curr_oprnd:
                case Operator():
                    # Uncast from tuple for list assembly
                    evald_oprnds.append(*curr_oprnd(args))
                case Expression():
                    evald_oprnds.append(curr_oprnd.val)
                case Argument():
                    # Python dict objects will raise a KeyError if no such key exists
                    evald_oprnds.append(args[curr_oprnd])
                case _:
                    raise TypeError
        # Cast/recast to tuple for data safety
        evald_oprnds = tuple(evald_oprnds)

        return self.evaluate(evald_oprnds)

    @abstractmethod
    def evaluate(
        self, 
        evald_oprnds : Tuple[SymbolicValueType, ...]
        ) -> Tuple[SymbolicValueType]:

        pass

class BinaryOperator(Operator, ABC):

    def __init__(
        self,
        first_oprnd : OperandInputType,
        second_oprnd : OperandInputType
        ):

        super().__init__([first_oprnd, second_oprnd])
    
    @abstractmethod
    def evaluate(
        self, 
        evald_oprnds : Tuple[SymbolicValueType, SymbolicValueType]
        ) -> Tuple[SymbolicValueType]:

        pass

class UnaryOperator(Operator, ABC):
    
    def __init__(
        self,
        oprnd : OperandInputType
        ):

        super().__init__([oprnd])
    
    @abstractmethod
    def evaluate(
        self,
        evald_oprnd : Tuple[SymbolicValueType]
        ) -> Tuple[SymbolicValueType]:

        pass

class SpatialOperator(UnaryOperator, ABC):
    """
    Spatial operators are best thought of as unary since they act in a single space.
    Even if their operand is vector-valued (or even tensor-valued), it is most accurately understood as a single entity, in a single space.
    As a result, there are no fundamental spatial operators which are non-unary.
    """

    def __init__(
        self,
        oprnd : OperandInputType,
        spce : Space  
        ):

        super().__init__(oprnd)
        self.spce = spce

    @abstractmethod
    def evaluate(
        self, 
        evald_oprnd : Tuple[SymbolicValueType]
        ) -> Tuple[SymbolicValueType]:
        
        pass


'''
Binary Operators
'''

class Sum(BinaryOperator):    

    def evaluate(
        self, 
        evald_oprnds : Tuple[SymbolicValueType, SymbolicValueType]
        ) -> Tuple[SymbolicValueType]:

        result = evald_oprnds[0] + evald_oprnds[1]
        return (result,)

class Difference(BinaryOperator):

    def evaluate(
        self, 
        evald_oprnds : Tuple[SymbolicValueType, SymbolicValueType]
        ) -> Tuple[SymbolicValueType]:

        result = evald_oprnds[0] - evald_oprnds[1]
        return (result,)

class ScalarMultiplication(BinaryOperator):

    def evaluate(
        self, 
        evald_oprnds : Union[
            Tuple[SymbolicScalarValueType, SymbolicValueType], 
            Tuple[SymbolicValueType, SymbolicScalarValueType]
            ]
        ) -> Tuple[SymbolicValueType]:

        result = evald_oprnds[0] * evald_oprnds[1]
        return (result,)

# Technically redundant due to ScalarMultiplication, but VERY convenient
class ScalarDivision(BinaryOperator):
    
    def evaluate(
        self, 
        evald_oprnds : Union[
            Tuple[SymbolicScalarValueType, SymbolicValueType], 
            Tuple[SymbolicValueType, SymbolicScalarValueType]
            ]
        ) -> Tuple[SymbolicValueType]:

        result = evald_oprnds[0] / evald_oprnds[1]
        return (result,)

class DotProduct(BinaryOperator):

    def evaluate(
        self, 
        evald_oprnds : Tuple[SymbolicVectorValueType, SymbolicVectorValueType]
        ) -> Tuple[SymbolicScalarValueType]:

        # SymPy Matrix objects have a .dot() method
        result = evald_oprnds[0].dot(evald_oprnds[1])
        return (result,)

class MatrixProduct(BinaryOperator):

    def evaluate(
        self,
        evald_oprnds : Tuple[SymbolicMatrixValueType, SymbolicMatrixValueType]
        ) -> Tuple[SymbolicMatrixValueType]:

        result = evald_oprnds[0] @ evald_oprnds[1]
        return (result,)


'''
Unary Operators
'''

class Negation(UnaryOperator):
    
    def evaluate(
        self,
        evald_oprnd : Tuple[SymbolicValueType]    
        ) -> Tuple[SymbolicValueType]:

        result = -evald_oprnd[0]
        
        return (result,)

class Derivative(UnaryOperator):

    def __init__(
        self,
        oprnd : OperandInputType,
        syms : List[sp.Symbol]
        ):

        super().__init__(oprnd)
        self.syms = syms

    def evaluate(
        self,
        evald_oprnd : Tuple[SymbolicValueType]
        ) -> Tuple[SymbolicValueType]:

        result = evald_oprnd[0]
        for curr_sym in self.syms:
            result = sp.Derivative(result, curr_sym)
        # Perform differentiation
        result = result.doit()

        return (result,)


'''
Spatial(Unary) Operators
'''

class Gradient(SpatialOperator):
    
    def evaluate(
        self,
        evald_oprnd : Tuple[SymbolicScalarValueType]
        ) -> Tuple[SymbolicVectorValueType]:

        result = sp.zeros(len(self.spce), 1)
        for curr_dim_idx in range(len(self.spce)):
            result[curr_dim_idx, 0] = sp.Derivative(evald_oprnd[0], self.spce[curr_dim_idx].sym)
        result = result.doit()

        return (result,)

class Divergence(SpatialOperator):

    def evaluate(
        self,
        evald_oprnd : Tuple[SymbolicVectorValueType]
        ) -> Tuple[SymbolicScalarValueType]:

        result = 0
        for curr_dim_idx in range(len(self.spce.dims)):
            result += sp.Derivative(evald_oprnd[0][curr_dim_idx], self.spce[curr_dim_idx].sym)
        result = result.doit()

        return (result,)

# For convenience
class Laplacian(SpatialOperator):

    def __init__(
        self,
        oprnd : OperandInputType,
        spce : Space 
        ):

        super().__init__(oprnd, spce)
        self.grad = Gradient(oprnd, spce)
        self.div = Divergence(oprnd, spce)

    def evaluate(
        self,
        evald_oprnd : Tuple[SymbolicScalarValueType]  
        ) -> Tuple[SymbolicScalarValueType]:
        
        return self.div.evaluate(self.grad.evaluate(evald_oprnd))

