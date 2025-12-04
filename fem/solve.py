# Libraries
import numpy as np


def solve(self, n_quad_points=2):

    glbl_num_vars = len(self.mesh.nds_vars_is)
    glbl_num_phys_els = len(self.mesh.els_nds_is)
    glbl_els_nds_is = self.mesh.els_nds_is
    glbl_nds_vec_crds = self.mesh.nds_vec_crds
    glbl_nds_vars_is = self.mesh.nds_vars_is

    # FEM matrices/vectors
    glbl_srcs = np.zeros([glbl_num_vars])
    glbl_op_coefs = np.zeros((glbl_num_vars, glbl_num_vars))

    for curr_el_i in range(glbl_num_phys_els):
        curr_el_nds_glbl_is = glbl_els_nds_is[curr_el_i]
        curr_el_nds_vars_glbl_is = np.concatenate([
            glbl_nds_vars_is[curr_nd_glbl_i] 
            for curr_nd_glbl_i in curr_el_nds_glbl_is
            ])
        curr_el_nds_vec_crds = glbl_nds_vec_crds[curr_el_nds_glbl_is]

        # --- Element contributions ---
        curr_el_vol, curr_el_bdry, curr_el_srcs = compute_element_contributions(
            (self.vol_funcs, self.bdry_funcs, self.src_funcs),
            curr_el_nds_vec_crds,
            n_points=n_quad_points
        )

        # Sum volume + boundary contributions
        curr_el_op_coefs = curr_el_vol + curr_el_bdry

        # --- Assemble into global matrices ---
        glbl_op_coefs[np.ix_(curr_el_nds_vars_glbl_is, curr_el_nds_vars_glbl_is)] += curr_el_op_coefs
        glbl_srcs[curr_el_nds_vars_glbl_is] += curr_el_srcs.flatten()

    soln = np.linalg.solve(glbl_op_coefs, glbl_srcs)
    soln = soln.flatten()