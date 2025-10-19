# Libraries
import numpy as np
from abc import ABC, abstractmethod


class BoundaryCondition:
    """
    General (base) boundary condition class.
    """
    types_registry = {}

    def __init_subclass__(cls, type:str=None, **kwargs):
        super().__init_subclass__(**kwargs)
        if type is not None:
            if type in BoundaryCondition.types_registry:
                raise ValueError(f"Duplicate boundary condition type: {type}")
            BoundaryCondition.types_registry[type] = cls

    def __init__(self, nds_glbl_is:np.ndarray, nds_vars_coefs:np.ndarray):
        self.nds_glbl_is = nds_glbl_is
        self.nds_vars_coefs = nds_vars_coefs
    
    @classmethod
    def create(cls, type:str, nds_glbl_is:np.ndarray, nds_vars_coefs:np.ndarray):
        """
        ## Purpose:
        User-facing method for creation of boundary conditions.

        ## Arguments:
        *Boundary condition types **must** be spelled correctly.*
        
        All nodes specified by `nds_glbl_is` must have the same stored variables (ex. nodes may not be a mix of those that store temperature T and pressure P). Define separate boundary condition instances in cases where you have groups of nodes with mismatched variables.
        | Name             |     | Data Type    |     | Description                                                                                                                                     |
        |:-----------------|-----|:-------------|-----|:------------------------------------------------------------------------------------------------------------------------------------------------|
        | `type`           |.....| `str`        |.....| Boundary condition type (ex. "Dirichlet")                                                                                                       |
        | `nds_glbl_is`    |.....| `np.ndarray` |.....| Vector of node-specifying indices associated with this boundary condition                                                                       |
        | `nds_vars_coefs` |.....| `np.ndarray` |.....| Array of coefficients specific to boundary condition type, used to calculate nodal variable values at that boundary. Shape: (num_nds, num_vars_per_nd) |

        ## Returns:
        A BoundaryCondition subclass instance (it's complicated).
        """

        if type not in cls.types_registry:
            raise ValueError(f"Unknown boundary condition type: {type}")
        return cls.types_registry[type](nds_glbl_is, nds_vars_coefs)

    @abstractmethod
    def apply(self, glbl_nds_vars_is:list, glbl_op_coefs:np.ndarray, glbl_srcs:np.ndarray) -> None:
        """
        ## Purpose:
        Applies a boundary condition to the FEM system.

        ## Arguments:
        Enforced arguments for boundary condition application functions are:
        | Name               |     | Data Type    |     | Description  |
        |:-------------------|-----|:-------------|-----|:-------------|
        | `glbl_nds_vars_is` |.....| `list`       |.....|              |
        | `glbl_op_coefs`    |.....| `np.ndarray` |.....|              |
        | `glbl_srcs`        |.....| `np.ndarray` |.....|              |

        ## Returns:
        Nothing.
        """

        pass

    def for_each_node(self, func, glbl_nds_vars_is):
        for curr_nd_glbl_i, curr_nd_vars_coefs in zip(self.nds_glbl_is, self.nds_vars_coefs):
            curr_nd_vars_glbl_is = glbl_nds_vars_is[curr_nd_glbl_i]
            func(curr_nd_vars_glbl_is, curr_nd_vars_coefs)

#-----
# Built-in boundary conditions
#-----


# Need to update to allow non-constant coefficients (probably supplied as sympy-created function handles?)
# apply() will need to allow function handles as input and will also need to take nodal coordinates as input
# also need to add support for non-hat weighting functions ):
class Dirichlet(BoundaryCondition, type="Dirichlet"):  
    """
    u = const.
    vars_coefs = [[c0], [c1], [c2]] to apply to all nodes with variablees [[v0], [v1], ...]
    """
    def apply(self, glbl_nds_vars_is:list, glbl_op_coefs:np.ndarray, glbl_srcs:np.ndarray):
        def apply_at_node(vars_glbl_is, vars_coefs):
            
            for curr_var_glbl_i, curr_var_coef in zip(vars_glbl_is, vars_coefs):
                glbl_op_coefs[curr_var_glbl_i, :              ] = 0.0 # hat weighting fn
                glbl_op_coefs[curr_var_glbl_i, curr_var_glbl_i] = 1.0 # hat weighting fn
                glbl_srcs[curr_var_glbl_i] = curr_var_coef
        self.for_each_node(apply_at_node, glbl_nds_vars_is)

class Neumann(BoundaryCondition, type="Neumann"):
    """
    u_n = const.
    vars_coefs = [[c0], [c1], [c2]] to apply to all nodes with variablees [[v0], [v1], ...]
    """
    def apply(self, glbl_nds_vars_is:list, glbl_op_coefs:np.ndarray, glbl_srcs:np.ndarray):
        def apply_at_node(vars_glbl_is, vars_coefs):
            for curr_var_glbl_i, curr_var_coef in zip(vars_glbl_is, vars_coefs):
                # assume hat weighting fns.
                glbl_srcs[curr_var_glbl_i] += curr_var_coef
        self.for_each_node(apply_at_node, glbl_nds_vars_is)

class Robin(BoundaryCondition, type="Robin"):
    """
    alpha*u + beta*u_n = gamma
    beta > 0, otherwise use Dirichlet BC
    vars_coefs = [[[alpha0], [beta0], [gamma0]], [[alpha1], [beta1], [gamma1]], ...] for vars = [[v0], [v1], ...]
    """
    def apply(self, glbl_nds_vars_is:list, glbl_op_coefs:np.ndarray, glbl_srcs:np.ndarray):
        def apply_at_node(vars_glbl_is, vars_coefs):
            for (curr_var_glbl_i, curr_var_coefs) in zip(vars_glbl_is, vars_coefs):
                alpha, beta, gamma = curr_var_coefs
                glbl_op_coefs[curr_var_glbl_i, curr_var_glbl_i] += (alpha / beta)
                glbl_srcs[curr_var_glbl_i] += (gamma / beta)
        self.for_each_node(apply_at_node, glbl_nds_vars_is)
