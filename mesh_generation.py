# Libraries
import numpy as np


'''
Intended workflow:

 user defines geometric regions
  |
  +-> mesh generator creates coupled (but independently stored and managed) connectivity matrices & nodal coordinates
       |
       +-> FEM uses both and adapts as needed

 - Nodal coordinates should be separated from connectivity to:
    - Reuse the same coordinates with different element types
    - Support adaptive meshing
    - Define boundaries using nodal coordinates
    - Support surface-specific curved elements

Tasks for curved boundaries/surfaces
1. Support for 1D, 2D, and 3D linear elements
2. Support for quadratic elements
3. Curved element edges only on curved boundaries, interior should stay linear
4. Functionality to project nodes exactly onto the boundary (spline?)
5. Optimization to precompute curved element characteristics and save runtime
'''

# uniform node spacing -- a lot of overlap with mesh_rectangle(), implies there should be object-oriented class structure here like with boundary conditions
def mesh_line(bdry_nds_vec_crds:np.ndarray, nums_nds):
    num_nds_x, = nums_nds

    # Use minimum point as starting point for consistent mesh orientation
    min_bdry_nd_i = np.argmin([np.sum(curr_nd_vec_crds) for curr_nd_vec_crds in bdry_nds_vec_crds])
    bdry_nds_vec_crds_ordered = np.concatenate((bdry_nds_vec_crds[min_bdry_nd_i:], bdry_nds_vec_crds[:min_bdry_nd_i]), axis=0)

    # Line vector definition
    x_vec = bdry_nds_vec_crds_ordered[1] - bdry_nds_vec_crds_ordered[0]

    # Nodal coordinate matrix & element-node connectivity matrix initialization
    nds_vec_crds = []
    els_nds_is = []

    # 1D map from nodal coordinate grid position (nd_x_i) to end-result node index
    nds_is = np.full((num_nds_x), fill_value=-1, dtype=int)
    num_nds_gend = 0

    # Nodal coordinate generation performed via linear interpolation
    for curr_nd_x_i in range(num_nds_x):
        curr_nd_x_vec = (curr_nd_x_i / (num_nds_x-1)) * x_vec

        # Nodal coordinate
        curr_nd_vec_crd = bdry_nds_vec_crds_ordered[0] + curr_nd_x_vec
        nds_vec_crds.append(curr_nd_vec_crd)

        # Nodal indexing
        nds_is[curr_nd_x_i] = num_nds_gend
        num_nds_gend += 1

        # Element specification
        if curr_nd_x_i > 0:
            curr_nd0_i = nds_is[curr_nd_x_i-1]
            curr_nd1_i = nds_is[curr_nd_x_i  ]
            # Line elements
            els_nds_is.append([curr_nd0_i, curr_nd1_i])

    nds_vec_crds = np.array(nds_vec_crds)
    els_nds_is = np.array(els_nds_is)
    return nds_vec_crds, els_nds_is

# Planar rectangle with orthogonal (90*) edges
def mesh_rectangle(bdry_nds_vec_crds:np.ndarray, nums_nds): # points are arranged in CCW order when looking down on rectangle in its plane
    num_nds_x, num_nds_y = nums_nds

    # Use minimum point as starting point for consistent mesh orientation
    min_bdry_nd_i = np.argmin([np.sum(curr_corner_crds) for curr_corner_crds in bdry_nds_vec_crds])
    bdry_nds_vec_crds_ordered = np.concatenate((bdry_nds_vec_crds[min_bdry_nd_i:], bdry_nds_vec_crds[:min_bdry_nd_i]), axis=0)

    # Rectangle edge vector definition
    edge_x_vec = bdry_nds_vec_crds_ordered[ 1] - bdry_nds_vec_crds_ordered[0]
    edge_y_vec = bdry_nds_vec_crds_ordered[-1] - bdry_nds_vec_crds_ordered[0]
    
    # Nodal coordinate matrix & element-node connectivity matrix initialization
    nds_vec_crds = []
    els_nds_is = []

    # 2D map from nodal coordinate grid position (nd_x_i, nd_y_i) to end-result node index
    nds_is = np.full((num_nds_x, num_nds_y), fill_value=-1, dtype=int)
    num_nds_gend = 0

    # Nodal coordinate grid generation performed via linear interpolation along rectangle edges
    for curr_nd_x_i in range(num_nds_x):
        curr_nd_x_vec = (curr_nd_x_i / (num_nds_x-1)) * edge_x_vec
        for curr_nd_y_i in range(num_nds_y):
            curr_nd_y_vec = (curr_nd_y_i / (num_nds_y-1)) * edge_y_vec

            # Nodal coordinate
            curr_nd_vec_crd = bdry_nds_vec_crds_ordered[0] + curr_nd_x_vec + curr_nd_y_vec
            nds_vec_crds.append(curr_nd_vec_crd)

            # Nodal indexing
            nds_is[curr_nd_x_i, curr_nd_y_i] = num_nds_gend
            num_nds_gend += 1

            # Element specification
            if curr_nd_x_i > 0 and curr_nd_y_i > 0: # Can't make any elements when there aren't enough nodes to do so
                curr_nd0_i = nds_is[curr_nd_x_i-1, curr_nd_y_i-1]
                curr_nd1_i = nds_is[curr_nd_x_i  , curr_nd_y_i-1]
                curr_nd2_i = nds_is[curr_nd_x_i  , curr_nd_y_i  ]
                curr_nd3_i = nds_is[curr_nd_x_i-1, curr_nd_y_i  ]
                # Rectangular elements
                els_nds_is.append([curr_nd0_i, curr_nd1_i, curr_nd2_i, curr_nd3_i])
    
    nds_vec_crds = np.array(nds_vec_crds)
    els_nds_is = np.array(els_nds_is)
    return nds_vec_crds, els_nds_is

# returns a list of np.ndarrays, to avoid wasted space (numpy arrays must be regular)
def gen_nds_vars_is(nums_vars_per_nd:np.ndarray, num_nds):
    nds_vars_is = [nums_vars_per_nd[curr_nd_i] * np.arange(curr_nd_i, curr_nd_i+1) for curr_nd_i in range(num_nds)]
    return nds_vars_is
