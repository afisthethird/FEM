# Libraries
import numpy as np
import matplotlib.pyplot as plt
# Scripts
from PYTHON.element_computes import compute_element_volume                   as comp_el_vol
from PYTHON.element_computes import compute_element_operator_coefficients    as comp_el_op_coefs
from PYTHON.element_computes import compute_element_sources                  as comp_el_srcs
from PYTHON.boundary_conditions import BoundaryCondition as BC

'''
Abbreviations Used:
    coef  : coefficient
    num   : number
    el    : element
    nd    : node
    var   : variable
    crd   : coordinate
    i     : index
    is    : indices
    val   : value
    bc/BC : boundary condition
    curr  : current
    src   : source
    op    : operator
    vol   : volume
    dim   : dimension
    DE    : differential equation
    deriv : derivative
    glbl  : global
    lcl   : local
    comp  : compute
    gen   : generate
    usr   : user
'''

'''
Variable Descriptions:

    num_els          : Number of elements in the system
    num_nds_per_el   : Number of nodes per element
    num_vars_per_nd  : Number of variables stored per node
        - Ex: storing temperature T and pressure P --> num_node_vars = 2
        - Typically referred to as "degrees of freedom" in FEM codes
    num_nds          : Number of nodes in the system
    num_vars         : Total number of stored variables across all system nodes

    glbl_nds_crds    : System nodal coordinate values
    glbl_els_nds_is  : System-wide element-node connectivity matrix
        - Each row specifies one element by index
        - Each column specifies that element's nodes by index
            - The first column specifies the start (min.) node
            - The second column specifies the end (max.) node
    
    bc_nds_is        : Specifies by index the nodes with imposed boundary conditions
    bc_nds_vars_vals : Specifies the imposed variable values at boundary nodes
        - Each column holds the imposed nodal variable values for the bc_nodes node at the same colum index

    glbl_srcs        : System-level vector containing discretized source terms
        - Represents external contributions (loads, heat generation, charge distribution, etc.) driving the system
        - Typically referred to as the "force vector" in FEM codes
    glbl_op_coefs    : System-level matrix containing discretized DE operator coefficients
        - Relates nodal variables (unkowns) to source terms in the system
        - Typically referred to as the "stiffness matrix" in FEM codes
'''


# returns a list of np.ndarrays, to avoid wasted space (numpy arrays must be regular)
def gen_nds_vars_is(nums_vars_per_nd:np.ndarray, num_nds):
    nds_vars_is = [nums_vars_per_nd[curr_nd_i] * np.arange(curr_nd_i, curr_nd_i+1) for curr_nd_i in range(num_nds)]
    return nds_vars_is

# Currently only supports line meshes
def gen_els_nds_is(nums_els:tuple) -> np.ndarray:
    """
    ## Purpose:
    Generates the element-to-node connectivity matrix for a specified element arrangement.
    
    ## Arguments:
    | Name       |     | Data Type |     | Description                                   |
    |:-----------|-----|:----------|-----|:----------------------------------------------|
    | `nums_els` |.....| `tuple`   |.....| Number of elements along each coordinate axis |

    ## Returns:
    | Name         |     | Data Type    |     | Description                                                         |
    |:-------------|-----|--------------|-----|:--------------------------------------------------------------------|
    | `els_nds_is` |.....| `np.ndarray` |.....| Element-node connectivity matrix, with shape (`num_els`, `num_nds`) |
    """

    num_dims = len(nums_els)
    
    if num_dims == 1:
        num_els = nums_els[0]
        els_nds_is = np.zeros((num_els, 2), dtype=int)
        els_nds_is[:, 0] = np.arange(0, num_els)
        els_nds_is[:, 1] = np.arange(1, num_els + 1)
    elif num_dims == 2:
        raise NotImplementedError("2D triangle mesh generation not implemented yet.")
    elif num_dims == 3:
        raise NotImplementedError("3D tetrahedron mesh generation not implemented yet.")
    else:
        raise ValueError(f"Unsupported number of dimensions: {num_dims}")
    
    return els_nds_is


# PDE specification
usr_a_coef =  1
usr_b_coef = -3
usr_c_coef =  2

# User-specified simulation parameters & global connectivity/indexing mesh generation
usr_nums_els = (3,)
glbl_els_nds_is = gen_els_nds_is(usr_nums_els)
glbl_num_els, glbl_num_nds_per_el = np.shape(glbl_els_nds_is) # NEEDS TO BE GENERALIZED!! I'm thinking a class for index meshes would be in order, with subclassed methods for 1D/2D/3D since behavior differs
glbl_num_nds = glbl_num_els + 1 # disgusting
usr_nums_vars_per_nd = np.ones(glbl_num_nds, dtype=int) # Need to generalize to custom variable count per node based on simulation regions; user should supply simulation regions that (auto) generate connectivity meshes
glbl_nds_vars_is = gen_nds_vars_is(usr_nums_vars_per_nd, glbl_num_nds)
glbl_num_vars = sum([curr_glbl_nd_vars_is.size for curr_glbl_nd_vars_is in glbl_nds_vars_is])

# Index generation
glbl_nds_crds = np.linspace(0, 1, glbl_num_els+1) # 1D, uniform spacing, x=[0,1] --- NEEDS to be generalized and combined with mesh generation

# Boundary conditions --- need to automate instead of manually specifying index numbers, again, should be by user-specified region
glbl_bcs = {
    "left":  BC.create("Dirichlet", np.array([0])             , np.array([[0.0]])),
    "right": BC.create("Dirichlet", np.array([glbl_num_nds-1]), np.array([[0.0]]))
}

# FEM matrices/vectors
glbl_srcs = np.zeros([glbl_num_vars])                    # [F]
glbl_op_coefs = np.zeros((glbl_num_vars, glbl_num_vars)) # [K]

for curr_el_i in range(glbl_num_els): # not currently suited to anything other than 1D, need to lock down stucture of element-node connectivity matrix (is it a list of numpy row vectors? unsure.)
    curr_el_nds_glbl_is = glbl_els_nds_is[curr_el_i]
    curr_el_nds_vars_glbl_is = np.concatenate([glbl_nds_vars_is[curr_nd_glbl_i] for curr_nd_glbl_i in curr_el_nds_glbl_is])

    curr_el_nds_crds = glbl_nds_crds[curr_el_nds_glbl_is].reshape(-1, 1)
    curr_el_vol = comp_el_vol(curr_el_nds_crds)

    curr_el_op_coefs = comp_el_op_coefs(usr_a_coef, usr_b_coef, usr_c_coef, curr_el_vol)
    curr_el_srcs = comp_el_srcs(curr_el_vol)

    glbl_op_coefs[np.ix_(curr_el_nds_vars_glbl_is, curr_el_nds_vars_glbl_is)] += curr_el_op_coefs
    glbl_srcs[curr_el_nds_vars_glbl_is] += curr_el_srcs.flatten()

for curr_bc_region, curr_bc in glbl_bcs.items():
    curr_bc.apply(glbl_nds_vars_is, glbl_op_coefs, glbl_srcs)
soln = np.linalg.solve(glbl_op_coefs, glbl_srcs)
soln = soln.flatten()




# Analytical solution coefficients
c1 = 0.5 / np.exp(1)
c2 = -0.5 * (1 + 1 / np.exp(1))

# Compute analytical solution at FEM nodes
xcor = glbl_nds_crds.copy()
esol = c1 * np.exp(2 * xcor) + c2 * np.exp(xcor) + 0.5

# Compute analytical solution at finer points for smooth plotting
x_fine = np.linspace(0, 1, 100)
esol_fine = c1 * np.exp(2 * x_fine) + c2 * np.exp(x_fine) + 0.5

# Print numerical and analytical solution at node locations
print(" Node   FEM_Solution   Exact_Solution")
for i in range(len(xcor)):
    print(f"{i+1:5d}   {soln[i]:.6f}       {esol[i]:.6f}")

# Plot FEM vs exact solution
plt.plot(xcor, soln, 'r-', linewidth=2, label='FEM')
plt.plot(xcor, esol, 'b-', linewidth=2, label='Exact (nodes)')
plt.plot(x_fine, esol_fine, 'k--', linewidth=2, label='Exact (fine)')
plt.xlabel('x')
plt.ylabel('u(x)')
plt.legend()
plt.grid(True)
plt.title('FEM vs Analytical Solution')
plt.show()