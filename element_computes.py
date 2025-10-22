# Libraries
import numpy as np
# Scripts
import PYTHON.reference_elements as ref_els

# drafting
def element_operator_integrand():

    # Compute operator terms for (i,j)
    diffusion = a_coef * 
    convection = b_coef * 
    reaction = c_coef * 

    return (-diffusion + convection + reaction) * 


# Currently limited to 1D DEs of the form: au"+bu'+cu = f(x) ; u = u(x)
# Currently limited to 1D hat functions
# messy at the moment
def compute_element_operator_coefficients(a_coef:float, b_coef:float, c_coef:float, el_nds_crds:np.ndarray) -> np.ndarray:
    """
    ## Purpose:
    Computes the (local) discretized DE operator coefficients matrix for an element, based on its volume.

    ## Arguments:
    | Name     |     | Data Type |     | Description                |
    |:---------|-----|:----------|-----|:---------------------------|
    | `a_coef` |.....| `float`   |.....| DE operator u" coefficient |
    | `b_coef` |.....| `float`   |.....| DE operator u' coefficient |
    | `c_coef` |.....| `float`   |.....| DE operator u coefficient  |
    
    ## Returns:
    | Name          |     | Data Type    |     | Description                                                                          |
    |:--------------|-----|:-------------|-----|:-------------------------------------------------------------------------------------|
    | `el_op_coefs` |.....| `np.ndarray` |.....| Discretized DE operator coefficients matrix for the specified PDE and element volume |
    """

    degree = len(el_nds_crds) - 1
    el_num_nds = degree + 1

    # Initialize element operator coefficients matrix
    el_op_coefs = np.zeros((el_num_nds, el_num_nds))

    for el_op_coefs_i in range(el_num_nds):
        for el_op_coefs_j in range(el_num_nds):
            func = lambda xi: element_operator_integrand(xi, el_op_coefs_i, el_op_coefs_j, a_coef, b_coef, c_coef, el_nds_crds)
            val, n_points_used = ref_els.adaptive_Gaussian_quadrature(func, np.array([[-1], [1]]))
            el_op_coefs[el_op_coefs_i, el_op_coefs_j] = val
    
    return el_op_coefs

# drafting
def element_sources_integrand(source_func=lambda x: 1.0):


    return f_val * 

# Currently limited to f(x) = 1
# Currently limited to 1D hat functions
def compute_element_sources(el_nds_crds:np.ndarray, source_func=lambda x:1.0) -> np.ndarray:
    """
    ## Purpose:
    Computes the (local) source vector for an element, based on its volume.

    ## Arguments:
    | Name     |     | Data Type |     | Description    |
    |:---------|-----|:----------|-----|:---------------|
    | `el_vol` |.....| `float`   |.....| Element volume |
    
    ## Returns:
    | Name      |     | Data Type    |     | Description                                 |
    |:----------|-----|:-------------|-----|:--------------------------------------------|
    | `el_srcs` |.....| `np.ndarray` |.....| Local, element-specific source/force vector |
    """

    degree = len(el_nds_crds) - 1
    el_num_nds = degree + 1

    # Initialize element operator coefficients matrix
    el_srcs = np.zeros(el_num_nds)

    for el_srcs_i in range(el_num_nds):
        func = lambda xi: element_sources_integrand(xi, el_srcs_i, el_nds_crds, source_func)
        val, n_points_used = ref_els.adaptive_Gaussian_quadrature(func, np.array([[-1], [1]]))
        el_srcs[el_srcs_i] = val
    
    return el_srcs
