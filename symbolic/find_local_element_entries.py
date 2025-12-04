import sympy as sp

# Symbolic variables
x, y, x_e, x_e1, y_e, y_e1, dx, dy = sp.symbols('x y x_e x_e1 y_e y_e1 dx dy')

# Bilinear shape functions at each node for axis-aligned, order-1, rectangular FEM elements
phi_0 = (x_e1 - x  ) * (y_e1 - y  ) / (dx * dy)
phi_1 = (x    - x_e) * (y_e1 - y  ) / (dx * dy)
phi_2 = (x    - x_e) * (y    - y_e) / (dx * dy)
phi_3 = (x_e1 - x  ) * (y    - y_e) / (dx * dy)

phis = [phi_0, phi_1, phi_2, phi_3]

# Partial derivative function for convenience
def partial_derivatives(phi):
    dphi_dx = sp.diff(phi, x)
    dphi_dy = sp.diff(phi, y)
    return dphi_dx, dphi_dy

# Stiffness matrix entry calculation
K_matrix = sp.zeros(4, 4)
for m in range(4):
    for n in range(4):
        dphi_m_dx, dphi_m_dy = partial_derivatives(phis[m])
        dphi_n_dx, dphi_n_dy = partial_derivatives(phis[n])
        integrand = (dphi_m_dx * dphi_n_dx + dphi_m_dy * dphi_n_dy)
        integral = sp.integrate(sp.integrate(integrand, (x, x_e, x_e1)), (y, y_e, y_e1))
        K_matrix[m, n] = integral

# Source vector entry calculation
F_vector = sp.zeros(1, 4)
for n in range(4):
    integral = sp.integrate(sp.integrate(phis[n], (x, x_e, x_e1)), (y, y_e, y_e1))
    F_vector[0, n] = integral

# Substitution to replace element nodal coordinates with dx & dy
subs_dict = {
    x_e  : 1,
    x_e1 : 1 + dx,
    y_e  : 1,
    y_e1 : 1 + dy
    }
K_matrix = K_matrix.subs(subs_dict)
F_vector = F_vector.subs(subs_dict)

# Output!!
sp.pprint(sp.simplify(K_matrix))
sp.pprint(sp.simplify(F_vector))