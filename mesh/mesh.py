# Libraries
# Scripts
from Code.types import *


# Delaunay/Advancing-Front meshing

# Tets are bad for incompressibility due to the "inf-sup (Ladyzhenskaya–Babuška–Brezzi) stability condition"
# In incompressible problems (Stokes flow, nearly incompressible elasticity), you need:
#   velocity / displacement space (V_h)   &   pressure space (Q_h)
# to satisfy the inf-sup condition:
#   {wacky equation}
# which essentially says that:
#   - velocity/displacement fields must be "rich enough"
#   - to represent divergence-free constraints
#   - for each admissible pressure mode
# and linear tets do not pass this condition, since they use "P1" displacement/velocity & constant pressure (P0) if using "mixed formulation"

# A bilinear surface is any surface that can be written as:
#   x(ξ,η) = A + Bξ + Cη + Dξη
# which includes flat quads, saddles, and twisted surfaces, but NOT general piecewise-curved hexes

'''
p-adaptive FEM textbooks:
    - Szabó & Babuška, Finite Element Analysis (1991)
    - Ainsworth & Oden, A Posteriori Error Estimation in FEM (2000)
Adaptive FEM software:
    - deal.II, libMesh, NGSolve: implement mixed conforming/non-conforming p-adaptivity based on error indicators.
    - Concept: mark for p-refinement, then decide whether to refine neighbors or insert constraints.
Practical rule:
    - Use conforming if local mesh modification is cheap; otherwise, use non-conforming with hanging node constraints.
'''

class Mesh:

    def __init__(self):
        return None