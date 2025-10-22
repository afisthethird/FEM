# Libraries
import numpy as np
import sympy as sp


# Kept seperate from the ReferenceElement class for now because this could be useful for physical element analysis later
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

# Currently only computes Lagrangian shape functions for 1D (linear) elements -- want to introduce Hermite
# For quad/hex: tensor-product calculation needed
# For tri/tet: barycentric (area) coordinates needed
def create_element_shape_functions_symbolic(vec_crd_syms, nds_vec_crds_syms, el_dim:int, el_poly_order:int):
    el_shape_fns_sym = []
    for curr_poly_order_i in range(el_poly_order+1):
        curr_shape_fn_numerator   = 1
        curr_shape_fn_denominator = 1
        curr_nd_i_vec_crd_syms = nds_vec_crds_syms[curr_poly_order_i]
        for curr_poly_order_j in range(el_poly_order+1):
            if curr_poly_order_j != curr_poly_order_i:
                curr_nd_j_vec_crd_syms = nds_vec_crds_syms[curr_poly_order_j]
                # Assume 1D (linear) elements, hence the [0] accesses
                curr_shape_fn_numerator   *= (vec_crd_syms[0]           - curr_nd_j_vec_crd_syms[0])
                curr_shape_fn_denominator *= (curr_nd_i_vec_crd_syms[0] - curr_nd_j_vec_crd_syms[0])
        el_shape_fns_sym.append(sp.simplify(curr_shape_fn_numerator / curr_shape_fn_denominator))
    return el_shape_fns_sym

# Currently only supports 1D (linear) elements
# Uses Gauss-Lobatto points, for well-placed interior nodes in elements with higher-order shape functions
def compute_element_interior_nodes_vector_coordinates(el_num_intr_nds:int, el_bdry_nds_vec_crds:np.ndarray):
    # Legendre polynomial of degree equal to 1D element interior nodes
    Legendre_polynomial = np.polynomial.legendre.Legendre.basis(
        deg=el_num_intr_nds, 
        domain=el_bdry_nds_vec_crds.flatten()
        )
    Legendre_polynomial_roots = Legendre_polynomial.roots()
    
    el_intr_nds_vec_crds = np.array(Legendre_polynomial_roots).reshape(-1, 1)
    return el_intr_nds_vec_crds

class ReferenceElement:
    def __init__(self, dim:int, num_bdry_nds:int, poly_order:int):
        if dim == 1 and num_bdry_nds == 2:
            self.bdry_nds_vec_crds = np.array([[-1], [1]])
            num_intr_nds = (poly_order + 1) - num_bdry_nds
            self.intr_nds_vec_crds = compute_element_interior_nodes_vector_coordinates(num_intr_nds, self.bdry_nds_vec_crds)
            self.nds_vec_crds = np.vstack((
                np.atleast_2d(self.bdry_nds_vec_crds[0]), 
                self.intr_nds_vec_crds, 
                np.atleast_2d(self.bdry_nds_vec_crds[1])
                ))
        else:
            raise NotImplementedError("Only 1D currently supported")
        
        self.dim = dim
        self.poly_order = poly_order
        self.vec_crd_syms = create_vector_coordinate_symbols(self.dim, f"ref_x")
        
        self.shape_fns_sym = None
    
    def create_shape_functions_symbolic(self, func):
        nds_vec_crds_syms = create_nodes_vector_coordinates_symbols(
            len(self.nds_vec_crds), 
            self.dim, 
            f"ref_nds_x"
            )
        temp_shape_fns_sym = func(self.vec_crd_syms, nds_vec_crds_syms, self.dim, self.poly_order)
        subs_dict = {
            nds_vec_crds_syms[curr_nd_i, 0]: float(self.nds_vec_crds[curr_nd_i, 0])
            for curr_nd_i in range(len(self.nds_vec_crds))
        }
        shape_fns_sym = np.array([
            sp.simplify(curr_shape_fn_sym.subs(subs_dict))
            for curr_shape_fn_sym in temp_shape_fns_sym
        ])
        self.shape_fns_sym = shape_fns_sym


# Helper function to create unique nodal coordinate variables for simultaneous handling of multiple nodes symbolically
def create_nodes_vector_coordinates_symbols(num_nds:int, num_scal_crds_per_nd_vec_crd:int, base_name:str):
    nds_vec_crds_syms = [
        [
            sp.Symbol(f"{base_name}_nd_{curr_nd_i}_crd_{curr_scal_crd_i}") 
            for curr_scal_crd_i in range(num_scal_crds_per_nd_vec_crd)
            ] 
        for curr_nd_i in range(num_nds)
        ]
    
    return np.array(nds_vec_crds_syms)

def create_vector_coordinate_symbols(num_scal_crds:int, base_name:str):
    vec_crd_syms = [
        sp.Symbol(f"{base_name}_crd_{curr_scal_crd_i}") 
        for curr_scal_crd_i in range(num_scal_crds)
        ]
    return vec_crd_syms

def create_reference_element_vector_coordinate_transformation_function_symbolic(ref_el:ReferenceElement, num_phys_dims_per_phys_nd_vec_crd:int):
    # Symbolic coordinates for physical nodes of a 1D element
    num_nds_per_el = len(ref_el.nds_vec_crds)
    phys_nds_vec_crds_syms = create_nodes_vector_coordinates_symbols(num_nds_per_el, num_phys_dims_per_phys_nd_vec_crd, f"phys_nds_x")

    # Mapping x(x_ref)
    ref_el_vec_crd_transform_fn_sym = sp.Matrix([
        sum(ref_el.shape_fns_sym[curr_el_nd_i] * phys_nds_vec_crds_syms[curr_el_nd_i][curr_vec_crd_dim_i]
        for curr_el_nd_i in range(num_nds_per_el))
        for curr_vec_crd_dim_i in range(num_phys_dims_per_phys_nd_vec_crd)
        ])

    return phys_nds_vec_crds_syms, ref_el_vec_crd_transform_fn_sym

def create_reference_element_Jacobian_function_symbolic(ref_el_vec_crd_transform_fn_sym, ref_vec_crd_syms):
    ref_el_Jacobian_fn_sym = ref_el_vec_crd_transform_fn_sym.diff(*ref_vec_crd_syms)
    return ref_el_Jacobian_fn_sym
