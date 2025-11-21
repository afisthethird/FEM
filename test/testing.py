# Libraries
import sympy as sp
import numpy as np
# Scripts
from Code.types import *
from Code.symbolic.space import R2, R2_ref
from Code.symbolic.math import Expression, Argument
from Code.symbolic.math import Derivative, Gradient, Divergence, Laplacian
from Code.symbolic.geometry import BoundarySurface, Domain


x, y = R2.dims_syms()

x_min_bdry = BoundarySurface('x_min', sp.Ge(x, -2), R2, Expression( 0))
x_max_bdry = BoundarySurface('x_max', sp.Le(x,  2), R2, Expression(10))
y_min_bdry = BoundarySurface('y_min', sp.Ge(y, -2), R2, Expression( 0))
y_max_bdry = BoundarySurface('y_max', sp.Le(y,  2), R2, Expression( 0))

center_dom = Domain('center', [x_min_bdry, x_max_bdry, y_min_bdry, y_max_bdry], R2)

test_crd = np.array([0, 0])
print(center_dom.contains(test_crd))

a = Expression(1)
b = Expression(sp.Matrix([2, 2]))
c = Expression(sp.sin(x) * sp.sin(y))
u = Argument('u')

L = a*Laplacian(u, R2) + b.dot(Gradient(u, R2)) + c*u

u_sub = sp.Function('u')(x, y)
L_eval = L({u : u_sub})[0]
sp.pprint(L_eval)

