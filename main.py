# Libraries
import numpy as np
import matplotlib.pyplot as plt
# Scripts
from PYTHON.element_computes import compute_element_volume                   as comp_el_vol
from PYTHON.element_computes import compute_element_operator_coefficients    as comp_el_op_coefs
from PYTHON.element_computes import compute_element_sources                  as comp_el_srcs
from PYTHON.boundary_conditions import BoundaryCondition as BC
import PYTHON.mesh_generation as mesh_gen

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

# PDE specification
usr_a_coef =  1
usr_b_coef = -3
usr_c_coef =  2

# Random crap, WIP
glbl_nds_crds, glbl_els_nds_is = mesh_gen.mesh_line([(0.0,), (1.0,)], (4,))
usr_nums_vars_per_nd = np.array([[1] for _ in range(len(glbl_nds_crds))]) # hacky

glbl_num_nds = len(glbl_nds_crds)
glbl_num_els = len(glbl_els_nds_is)
glbl_nds_vars_is = mesh_gen.gen_nds_vars_is(usr_nums_vars_per_nd, glbl_num_nds)
glbl_num_vars = sum([curr_glbl_nd_vars_is.size for curr_glbl_nd_vars_is in glbl_nds_vars_is])

# Boundary conditions --- need to automate instead of manually specifying index numbers, again, should be by user-specified region
# Also need to update variable coefficients vector to a list (in line with custom numbers of variables per node)
glbl_bcs = {
    "left":  BC.create("Dirichlet", np.array([0])             , np.array([[0]])),
    "right": BC.create("Dirichlet", np.array([glbl_num_nds-1]), np.array([[0]]))
}

# FEM matrices/vectors
glbl_srcs = np.zeros([glbl_num_vars])
glbl_op_coefs = np.zeros((glbl_num_vars, glbl_num_vars))

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


# Plot FEM vs exact solution
plt.plot(glbl_nds_crds, soln, 'r-', linewidth=2, label='FEM')
plt.xlabel('x')
plt.ylabel('u(x)')
plt.legend()
plt.grid(True)
plt.title('FEM Solution')
plt.show()
