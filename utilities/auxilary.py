# Libraries
import sympy as sp
import numpy as np
from dataclasses import dataclass
# Scripts
from Code.types import *


# Utility/wrapper class for storing values associated with symbols
@dataclass
class StoredVariable():

    sym : sp.Symbol
    val : NumericValueType = None

    def __repr__(self):

        result = ""
        if self.val is not None:
            result += f"{self.val}"
        else:
            result += f"{self.sym}"

        return result


def make_callable(
    syms : Tuple[sp.Symbol, ...],
    expr : sp.Expr
    ) -> Callable[..., NumericValueType]:

    result = sp.lambdify(
        args = syms, 
        expr = expr, 
        modules = 'numpy'
        )

    return result


'''
Random crap from my previous architecture that I've kept around in case it's useful
'''

# Helper function to create unique coordinate variables for spaces
def construct_vector_coordinates_symbols(
    base_name   : str,
    num_of_dims : int,
    ) -> Tuple[sp.Symbol, ...]:
    
    return tuple(
        sp.Symbol(f"{base_name}_crd_{curr_dim_i}") 
        for curr_dim_i in range(num_of_dims)
        )


# Constructs nested lists of strings for coordinate-agnostic vertex/node indexing
def generate_vertices_indices_names(base_name, num_of_topo_dims, nums_of_verts_per_topo_dim):

    def recursively(curr_topo_dim, curr_vert_topo_dims_is):
        # Base case: Indices are generated as coordinate placeholders
        if curr_topo_dim == num_of_topo_dims:
            return (
                f"{base_name}" + 
                "".join([f"{curr_vert_topo_dim_i}" for curr_vert_topo_dim_i in curr_vert_topo_dims_is])
                )
        # Recursive case: Build a containing list for this topological dimension's axis
        return [
            recursively(curr_topo_dim+1, curr_vert_topo_dims_is+[curr_topo_dim_curr_num_of_verts])
                for curr_topo_dim_curr_num_of_verts in range(nums_of_verts_per_topo_dim[curr_topo_dim])
            ]

    return recursively(0, [])


def compute_Gauss_Lobatto_values(
    num_of_vals : int,
    start_val   : Union[int, float],
    end_val     : Union[int, float],
    ) -> np.ndarray:
    if num_of_vals == 2:
        return np.array([start_val, end_val], dtype=float)
    Legendre_poly = np.polynomial.legendre.Legendre.basis(deg=num_of_vals-1, domain=[start_val, end_val])
    Legendre_poly_roots = Legendre_poly.deriv().roots()
    vals = np.concatenate(([start_val], Legendre_poly_roots, [end_val]))
    return vals


# Keep around because possibly useful for analysis later
def compute_element_volume(el_dim:int, el_nds_crds:np.ndarray) -> float:
    num_el_nds = len(el_nds_crds)
    
    el_vol = 0

    if el_dim == 1 and num_el_nds == 2:
        # 1D line element (length)
        a, b = el_nds_crds
        el_vol = np.linalg.norm(b - a)

    elif el_dim == 2 and num_el_nds == 3:
        # 2D triangle element (area)
        a, b, c = el_nds_crds
        el_vol = 0.5 * np.linalg.norm(np.cross((b - a), (c - a)))

    elif el_dim == 3 and num_el_nds == 4:
        # 3D tetrahedron element (volume)
        a, b, c, d = el_nds_crds
        el_vol = np.abs(np.dot((a - d), np.cross((b - d), (c - d)))) / 6.0

    else:
        raise ValueError("Invalid element shape or dimensionality")
    
    return el_vol

