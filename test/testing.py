# Libraries
import sympy as sp
import numpy as np
# Scripts
from Code.types import *
from Code.symbolic.space import R2
from Code.symbolic.math import Expression, Argument
from Code.symbolic.math import Derivative, Gradient, Divergence, Laplacian
from Code.symbolic.geometry import Boundary, Domain
from Code.fem.equation import GoverningEquation


x, y = R2.dims_syms()

x_min_bdry = Boundary('x_min', sp.Ge(x, 0.0), R2)
x_max_bdry = Boundary('x_max', sp.Le(x, 1.0), R2)
y_min_bdry = Boundary('y_min', sp.Ge(y, 0.0), R2)
y_max_bdry = Boundary('y_max', sp.Le(y, 1.0), R2)

center_dom = Domain('center', [x_min_bdry, x_max_bdry, y_min_bdry, y_max_bdry], R2)

V = Argument('V')
ρ = 0.0
ɛ0 = 8.8541878188e-12 # vacuum
ɛr = 1.0 # Air
ɛ = ɛr/ɛ0

op = Laplacian(V, R2)
src = Expression(-ρ / ɛ)

eqn = GoverningEquation(
    name = 'Poisson',
    op = op,
    src = src,
    host_spce = R2
    )

w = Argument('w')
eqn.construct_weak_form_integrands(wghtng_plchldr = w)

V_sub = sp.Function('V')(x, y)
w_sub = sp.Function('w')(x, y)
sp.pprint(eqn.weak_form_op_intgrnd({V: V_sub, w: w_sub})[0])
sp.pprint(eqn.weak_form_src_intgrnd({w: w_sub})[0])