# Libraries
import numpy as np
from abc import ABC, abstractmethod


class BoundaryCondition:
    """"
    General (base) boundary condition class.
    """
    types_registry = {}

    def __init_subclass__(cls, type:str=None, **kwargs):
        super().__init_subclass__(**kwargs)
        if type is not None:
            if type in BoundaryCondition.types_registry:
                raise ValueError(f"Duplicate boundary condition type: {type}")
            BoundaryCondition.types_registry[type] = cls

    def __init__(self, nds_glbl_is:np.ndarray, nds_vars_vals:np.ndarray):
        self.nds_glbl_is = nds_glbl_is
        self.nds_vars_vals = nds_vars_vals
    
    @classmethod
    def create(cls, type:str, nds_is:np.ndarray, nds_vars_vals:np.ndarray):
        """
        ## Purpose:
        User-facing method for creation of boundary conditions.

        ## Arguments:
        *Boundary condition types **must** be spelled correctly.*
        | Name            |     | Data Type    |     | Description                                                                                |
        |:----------------|-----|:-------------|-----|:-------------------------------------------------------------------------------------------|
        | `type`          |.....| `str`        |.....| Boundary condition type (ex. "Dirichlet")                                                  |
        | `nds_is`        |.....| `np.ndarray` |.....| List of node-specifying indices associated with this boundary condition                    |
        | `nds_vars_vals` |.....| `np.ndarray` |.....| List of values for variables at each specified node, with shape (num_nds, num_vars_per_nd) |

        ## Returns:
        A BoundaryCondition subclass instance (it's complicated).
        """

        if type not in cls.types_registry:
            raise ValueError(f"Unknown boundary condition type: {type}")
        return cls.types_registry[type](nds_is, nds_vars_vals)

    @abstractmethod
    def apply(self, glbl_nds_vars_is:list, glbl_op_coefs:np.ndarray, glbl_srcs:np.ndarray) -> None:
        """
        ## Purpose:
        Applies a boundary condition to the FEM system.

        ## Arguments:
        Enforced arguments for boundary condition application functions are:
        | Name            |     | Data Type    |     | Description  |
        |:----------------|-----|:-------------|-----|:-------------|
        | `glbl_op_coefs` |.....| `np.ndarray` |.....|              |
        | `glbl_srcs`     |.....| `np.ndarray` |.....|              |

        ## Returns:
        Nothing.
        """

        pass

    def for_each_node(self, func, glbl_nds_vars_is):
        for curr_nd_glbl_i, curr_nd_vars_vals in zip(self.nds_glbl_is, self.nds_vars_vals):
            curr_nd_vars_glbl_is = glbl_nds_vars_is[curr_nd_glbl_i]
            func(curr_nd_vars_glbl_is, curr_nd_vars_vals)

#-----
# Built-in boundary conditions
#-----

class Dirichlet(BoundaryCondition, type="Dirichlet"):
    def apply(self, glbl_nds_vars_is:list, glbl_op_coefs:np.ndarray, glbl_srcs:np.ndarray):
        def apply_at_node(vars_glbl_is, vars_vals):
            for curr_var_glbl_i, curr_var_val in zip(vars_glbl_is, vars_vals):
                glbl_op_coefs[curr_var_glbl_i, :              ] = 0.0
                glbl_op_coefs[curr_var_glbl_i, curr_var_glbl_i] = 1.0
                glbl_srcs[curr_var_glbl_i] = curr_var_val
        self.for_each_node(apply_at_node, glbl_nds_vars_is)