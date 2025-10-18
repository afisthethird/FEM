# Libraries
import numpy as np


def compute_element_volume(el_crds:np.ndarray) -> float:
    """
    ## Purpose:
    Computes the space occupied by an element, defined by its nodal coordinates.
    - 1D & 2 nodes : line, returns length
    - 2D & 3 nodes : triangle, returns area
    - 3D & 4 nodes : tetrahedron, returns volume

    ## Arguments:
    *Coordinates **must** be supplied in Cartesian format.*
    | Name          |     | Data Type    |     | Description                                                                         |
    |:--------------|-----|:-------------|-----|:------------------------------------------------------------------------------------|
    | `el_nds_crds` |.....| `np.ndarray` |.....| `ndarray` of nodal coordinates for an element, with shape(`num_el_nds`, `num_dims`) |
    
    ## Returns:
    | Name     |     | Data Type |     | Description                                             |
    |:---------|-----|:----------|-----|:--------------------------------------------------------|
    | `el_vol` |.....| `float`   |.....| element volume, based on the supplied nodal coordinates |
    """
    
    num_nds_per_el, num_dims = el_crds.shape
    el_vol = None

    if num_dims == 1 and num_nds_per_el == 2:
        # 1D line element
        x1 = el_crds[0, 0]
        x2 = el_crds[1, 0]
        el_vol = np.abs(x2 - x1)
    elif num_dims == 2 and num_nds_per_el == 3:
        # 2D triangle element
        x1, y1 = el_crds[0]
        x2, y2 = el_crds[1]
        x3, y3 = el_crds[2]
        el_vol = 0.5 * np.abs((x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1))
    elif num_dims == 3 and num_nds_per_el == 4:
        # 3D tetrahedron element
        a, b, c, d = el_crds
        el_vol = np.abs(np.dot((a - d), np.cross((b - d), (c - d)))) / 6.0
    else:
        raise ValueError("Invalid element shape or dimensionality")
    
    return el_vol


# Currently limited to 1D DEs of the form: au"+bu'+cu = f(x) ; u = u(x)
# Currently limited to 1D hat functions
def compute_element_operator_coefficients(a_coef:float, b_coef:float, c_coef:float, el_vol:float) -> np.ndarray:
    """
    ## Purpose:
    Computes the (local) discretized DE operator coefficients matrix for an element, based on its volume.

    ## Arguments:
    | Name     |     | Data Type |     | Description                |
    |:---------|-----|:----------|-----|:---------------------------|
    | `a_coef` |.....| `float`   |.....| DE operator u" coefficient |
    | `b_coef` |.....| `float`   |.....| DE operator u' coefficient |
    | `c_coef` |.....| `float`   |.....| DE operator u coefficient  |
    | `el_vol` |.....| `float`   |.....| Element volume             |
    
    ## Returns:
    | Name          |     | Data Type    |     | Description                                                                          |
    |:--------------|-----|:-------------|-----|:-------------------------------------------------------------------------------------|
    | `el_op_coefs` |.....| `np.ndarray` |.....| Discretized DE operator coefficients matrix for the specified PDE and element volume |
    """

    diffusion_op_coefs  = (-a_coef / el_vol    ) * np.array([[ 1, -1], [-1,  1]])
    convection_op_coefs = ( b_coef          / 2) * np.array([[-1,  1], [-1,  1]])
    reaction_op_coefs   = ( c_coef * el_vol / 6) * np.array([[ 2,  1], [ 1,  2]])
    el_op_coefs = diffusion_op_coefs + convection_op_coefs + reaction_op_coefs
    return el_op_coefs


# Currently limited to f(x) = 1
# Currently limited to 1D hat functions
def compute_element_sources(el_vol:float) -> np.ndarray:
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

    el_srcs = np.array([(el_vol / 2), (el_vol / 2)])
    return el_srcs
