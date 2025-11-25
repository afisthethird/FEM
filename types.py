# Libraries
import sympy as sp
import numpy as np
from typing import TypeAlias, TypeVar, ClassVar, Any
from typing import Union, Annotated, Literal
from typing import Type, Self, Callable, List, Tuple, Dict
from typing import get_origin, get_args


# Utility function for extracting concrete, base type objects from aliased, non-base type objects
def extract_base_types(type_obj : Any) -> Tuple[type, ...]:

    # Extract nested type objects
    type_obj_origin = get_origin(type_obj)
    type_obj_args = get_args(type_obj)
    # Already a base type
    if type_obj_origin is None:
        # Avoid including Literals in extracted types
        if isinstance(type_obj, type):
            return (type_obj,)
        return ()
    # Annotation stripping
    elif type_obj_origin is Annotated:
        return extract_base_types(type_obj_args[0])
    # Non-base type
    else:
        result = ()
        for curr_type_obj_arg in type_obj_args:
            result += extract_base_types(curr_type_obj_arg)
        return result


# Shape setup
RowVectorValueShape    : TypeAlias = Literal["(m, 1)"]
ColumnVectorValueShape : TypeAlias = Literal["(1, n)"]
VectorValueShape       : TypeAlias = Union[
    RowVectorValueShape,
    ColumnVectorValueShape
    ]
MatrixValueShape       : TypeAlias = Literal["(m, n)"]
TensorValueShape       : TypeAlias = Literal["(m, n, ..)"]


# Centralized handling of typing regarding SymPy's convoluted matrix/tensor classes
SymPyMatrixType : TypeAlias = sp.MatrixBase
SymPyArrayType : TypeAlias = sp.NDimArray


# Supported symbolic value types
SymbolicScalarValueType : TypeAlias = Union[
    int, 
    float, 
    sp.Expr
    ]
SymbolicRowVectorValueType : TypeAlias = Annotated[
    SymPyMatrixType, 
    RowVectorValueShape
    ]
SymbolicColumnVectorValueType : TypeAlias = Annotated[
    SymPyMatrixType, 
    ColumnVectorValueShape
    ]
SymbolicVectorValueType : TypeAlias = Union[
    SymbolicRowVectorValueType, 
    SymbolicColumnVectorValueType
    ]
SymbolicMatrixValueType : TypeAlias = Annotated[
    SymPyMatrixType, 
    MatrixValueShape
    ]
SymbolicTensorValueType : TypeAlias = Annotated[
    SymPyArrayType,
    TensorValueShape
    ]
SymbolicValueType : TypeAlias = Union[
    SymbolicScalarValueType, 
    SymbolicVectorValueType, 
    SymbolicMatrixValueType,
    SymbolicTensorValueType
    ]


# Centralized handling of NumPy classes
NumPyIntegerType : TypeAlias = np.integer
NumPyDecimalType : TypeAlias = np.floating
NumPyScalarType : TypeAlias = Union[NumPyIntegerType, NumPyDecimalType]


# Supported numeric value types
NumericIntegerValueType : TypeAlias = Union[int, NumPyIntegerType]
NumericDecimalValueType : TypeAlias = Union[float, NumPyDecimalType]
NumericScalarValueType : TypeAlias = Union[
    NumericIntegerValueType,
    NumericDecimalValueType
    ]
NumericRowVectorValueType : TypeAlias = np.ndarray[
    NumericScalarValueType,
    RowVectorValueShape
    ]
NumericColumnVectorValueType : TypeAlias = np.ndarray[
    NumericScalarValueType,
    ColumnVectorValueShape
    ]
NumericVectorValueType : TypeAlias = Union[
    NumericRowVectorValueType, 
    NumericColumnVectorValueType
    ]
NumericMatrixValueType : TypeAlias = np.ndarray[
    NumericScalarValueType, 
    MatrixValueShape
    ]
NumericTensorValueType : TypeAlias = np.ndarray[
    NumericScalarValueType,
    TensorValueShape
    ]
NumericValueType : TypeAlias = Union[
    NumericScalarValueType, 
    NumericVectorValueType, 
    NumericMatrixValueType,
    NumericTensorValueType
    ]


# Combined types
ScalarValueType : TypeAlias = Union[
    SymbolicScalarValueType, 
    NumericScalarValueType
    ]
VectorValueType : TypeAlias = Union[
    SymbolicVectorValueType,
    NumericVectorValueType
    ]
MatrixValueType : TypeAlias = Union[
    SymbolicMatrixValueType,
    NumericMatrixValueType
    ]
TensorValueType : TypeAlias = Union[
    SymbolicTensorValueType,
    NumericTensorValueType
    ]
ValueType : TypeAlias = Union[
    ScalarValueType,
    VectorValueType,
    MatrixValueType,
    TensorValueType
    ]


# Non-numeric/symbolic common setup
IndexType : TypeAlias = NumericIntegerValueType
NameType: TypeAlias = str


# Runtime-usable extractions
# Runtime-usable extractions
RUNTIME_SYMBOLIC_SCALAR_VALUE_TYPE = extract_base_types(SymbolicScalarValueType)
RUNTIME_SYMBOLIC_VECTOR_VALUE_TYPE = extract_base_types(SymbolicVectorValueType)
RUNTIME_SYMBOLIC_MATRIX_VALUE_TYPE = extract_base_types(SymbolicMatrixValueType)
RUNTIME_SYMBOLIC_TENSOR_VALUE_TYPE = extract_base_types(SymbolicTensorValueType)
RUNTIME_INDEX_TYPE = extract_base_types(IndexType)
RUNTIME_NAME_TYPE = extract_base_types(NameType)

