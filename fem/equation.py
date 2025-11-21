# Libraries
import sympy as sp
import numpy as np
# Scripts
import PYTHON.symbolic as smblc
import PYTHON.element as elmnt
import PYTHON.space as spce
import PYTHON.auxilary as aux


class FiniteElementMethodEquation():
    def __init__(
        self,
        name                   : str,
        space                  : spce.Space,
        op                     : smblc.Operator,
        src_expr               : sp.Expr,
        el_poly_orders_per_dim : tuple[int]
        ):

        self.name = name

        self.phys_space = space
        self.ref_space = self.construct_reference_space()

        self.phys_op = op
        self.phys_src_expr = src_expr
        
        self.template_el = elmnt.TemplateElement(
            name                = 'el', 
            phys_space          = self.phys_space,
            ref_space           = self.ref_space,
            poly_orders_per_dim = el_poly_orders_per_dim
            )

        # Placeholders used for weak form conversion and coordinate transformation
        self.placeholder_test_fn_expr = sp.Symbol('w')
        self.placeholder_bdry_norm_vec_expr = sp.Matrix(sp.symbols(f'n0:{self.phys_space.num_of_dims}'))

        self.phys_weak_op_vol_part, self.phys_weak_op_bdry_part, self.phys_weak_src_expr = self.convert_to_weak_form()

        self.ref_weak_op_vol_part, self.ref_weak_op_bdry_part, self.ref_weak_src_expr = self.transform_physical_to_reference()

    def construct_reference_space(self):
        num_of_dims = self.phys_space.num_of_dims
        verts_vec_crds_vals = np.stack(
            np.meshgrid(
                *([[-1.0, 1.0]] * num_of_dims),
                indexing = 'ij'
                ),
            axis = -1
            )
        return spce.ReferenceSpace('ref_space', self.phys_space, verts_vec_crds_vals)

    # Need to avoid weak form conversion post-coordinate-transformation (computationally intense)
    # Optional perform_int_by_parts flag to support later to support later development of higher-order shape functions or storage of solution derivatives
    # Currently doesn't support general DifferentialOperators
    def convert_to_weak_form(self, perform_int_by_parts=True):  

        strong_op = self.phys_op
        strong_src_expr = self.phys_src_expr
        space = self.phys_space

        # Predefined differential Operators for easy integration by parts
        grad_op = smblc.Gradient(space)

        # Placeholder test function for weak form conversion
        placeholder_test_fn_expr = self.placeholder_test_fn_expr
        placeholder_test_fn_op = smblc.Expression(space, placeholder_test_fn_expr)
        # Placeholder boundary normal vector for integration by parts
        placeholder_bdry_norm_vec_expr = self.placeholder_bdry_norm_vec_expr

        # Recursive helper to perform abstract integration by parts
        # Only has support for the Laplacian at the moment, but it wouldn't be too hard to add support for general order-2 integration-by-parts or higher-order differential terms
        # ∫(∇²u(x)·v(x))dΩ = -∫(∇u·∇v)dΩ + ∫(∇u·vŜ)dS ; S = ∂Ω
        def integrate_by_parts(op:smblc.Operator):

            op_vol_part = None
            op_bdry_part = None

            # Examine components of CompositeOperators
            if isinstance(op, smblc.CompositeOperator):
                # Look for Laplacian: Divergence @ Gradient
                if op.type == 'of' and (isinstance(op.first_op, smblc.Divergence) and isinstance(op.second_op, smblc.Gradient)):
                    # Volume term: -∫(∇u·∇v)dΩ
                    op_vol_part = -1 * grad_op.of(placeholder_test_fn_op).dot(grad_op)
                    # Boundary term: ∫(∇u·vŜ)dS
                    op_bdry_part = grad_op.dot(placeholder_test_fn_op * placeholder_bdry_norm_vec_expr)
                # No Laplacian found
                else:
                    # Common to all composition types
                    first_op_vol_part, first_op_bdry_part = integrate_by_parts(op.first_op)
                    second_op_vol_part, second_op_bdry_part = integrate_by_parts(op.second_op)
                    op_vol_part = smblc.CompositeOperator(
                        space     = space,
                        first_op  = first_op_vol_part,
                        second_op = second_op_vol_part,
                        type      = op.type
                        )
                    op_bdry_part = smblc.CompositeOperator(
                        space     = space,
                        first_op  = first_op_bdry_part,
                        second_op = second_op_bdry_part,
                        type      = op.type
                        )
            # Non-Operator-involved Operators
            elif isinstance(op, (smblc.ScalarProduct, smblc.VectorProduct, smblc.MatrixProduct)):
                inner_op_vol_part, inner_op_bdry_part = integrate_by_parts(op.inner_op)
                coef_expr = op.expr
                if isinstance(op, smblc.ScalarProduct):
                    op_vol_part = coef_expr * inner_op_vol_part
                    op_bdry_part = coef_expr * inner_op_bdry_part
                elif isinstance(op, smblc.VectorProduct):
                    op_vol_part = inner_op_vol_part.dot(coef_expr)
                    op_bdry_part = inner_op_bdry_part.dot(coef_expr)
                elif isinstance(op, smblc.MatrixProduct):
                    op_vol_part = coef_expr @ inner_op_vol_part
                    op_bdry_part = coef_expr @ inner_op_bdry_part
            # Gradient needs special handling since it is vector-valued
            elif isinstance(op, smblc.Gradient):
                op_vol_part = op * placeholder_test_fn_op
                op_bdry_part = smblc.VectorZero(space)
            # Everything else can be left alone since it's assumed order-1 or less
            else:
                op_vol_part = op * placeholder_test_fn_op
                op_bdry_part = smblc.ScalarZero(space)
            
            return op_vol_part, op_bdry_part

        weak_op_vol_part, weak_op_bdry_part = integrate_by_parts(strong_op)
        weak_src_expr = strong_src_expr * placeholder_test_fn_expr

        return weak_op_vol_part, weak_op_bdry_part, weak_src_expr

    def transform_physical_to_reference(self):

        # Spaces
        phys_space = self.phys_space
        ref_space = self.ref_space

        # Placeholder test function
        ref_placeholder_test_fn_op = smblc.Expression(ref_space, self.placeholder_test_fn_expr)

        # Predefined Operators for easy substitution
        ref_grad_op = smblc.Gradient(ref_space)
        ref_div_op = smblc.Divergence(ref_space)
        ref_id_op = smblc.Identity(ref_space)

        # Vector coordinate symbols
        phys_vec_crd_syms = phys_space.vec_crd_syms

        # FEM Equation components
        phys_weak_op_vol_part = self.phys_weak_op_vol_part
        phys_weak_op_bdry_part = self.phys_weak_op_bdry_part
        phys_weak_src_expr = self.phys_weak_src_expr

        # Template element coordinate transformation expressions
        vec_crd_transform_expr = self.template_el.ref_to_phys_vec_crd_transform_expr
        Jacobian_expr = self.template_el.ref_to_phys_Jacobian_expr
        Jacobian_inv_expr = Jacobian_expr.inv()
        Jacobian_inv_transpose_expr = Jacobian_inv_expr.T
        Jacobian_det_expr = Jacobian_expr.det()

        # Substitution dictionary for all physical vector coordinate Symbols
        subs_dict = {
            phys_vec_crd_syms[curr_dim_i]: vec_crd_transform_expr[curr_dim_i]
            for curr_dim_i in range(self.phys_space.num_of_dims)
            }

        # Recursive helper for substituting into Expressions
        def transform_expression(phys_expr:sp.Expr):
            # Catch non-SymPy components
            if not isinstance(phys_expr, sp.Basic):
                return phys_expr
            # Substitute into any SymPy Expressions
            ref_expr = phys_expr.subs(subs_dict)
            # Only recurse into non-atomic Expressions (those with arguments)
            if ref_expr.args:
                # Transform arguments inside any Functions
                ref_expr = ref_expr.func(
                    *(transform_expression(curr_ref_arg)
                        for curr_ref_arg in ref_expr.args
                        )
                    )
            # Substitute return
            return ref_expr

        # Recursive helper for chain rule chaining
        def transform_operator(phys_op:smblc.Operator):
            # Components of CompositeOperators
            if isinstance(phys_op, smblc.CompositeOperator):
                return smblc.CompositeOperator(
                    space     = ref_space,
                    first_op  = transform_operator(phys_op.first_op), 
                    second_op = transform_operator(phys_op.second_op), 
                    type      = phys_op.type
                    )
            # Non-Operator-involved Operators
            elif isinstance(phys_op, (smblc.ScalarProduct, smblc.VectorProduct)):
                ref_inner_op = transform_operator(phys_op.inner_op)
                ref_coef_expr = transform_expression(phys_op.expr)
                if isinstance(phys_op, smblc.ScalarProduct):
                    return smblc.ScalarProduct(ref_space, ref_coef_expr, ref_inner_op)
                elif isinstance(phys_op, smblc.VectorProduct):
                    return smblc.VectorProduct(ref_space, ref_coef_expr, ref_inner_op)
                elif isinstance(phys_op, smblc.MatrixProduct):
                    return smblc.MatrixProduct(ref_space, ref_coef_expr, ref_inner_op)
            # Gradient: ∇f(x) = J⊥·∇f(ξ)
            if isinstance(phys_op, smblc.Gradient):
                return ref_grad_op @ Jacobian_inv_transpose_expr
            # Divergence: ∇·v(x) = (1/|J|)*∇·(|J|*J·v(ξ))
            elif isinstance(phys_op, smblc.Divergence):
                return (1 / Jacobian_det_expr) * ref_div_op.of(smblc.Identity @ (Jacobian_det_expr * Jacobian_expr))
            # Identity Operator
            elif isinstance(phys_op, smblc.Identity):
                return ref_id_op
            # Placeholder Operator
            elif isinstance(phys_op, smblc.Expression):
                return ref_placeholder_test_fn_op
            # Non-Operators
            else:
                return transform_expression(phys_op)
            
        ref_weak_op_vol_part = transform_operator(phys_weak_op_vol_part)
        ref_weak_op_bdry_part = transform_operator(phys_weak_op_bdry_part)
        ref_weak_src_expr = transform_expression(phys_weak_src_expr)

        return ref_weak_op_vol_part, ref_weak_op_bdry_part, ref_weak_src_expr

    def generate_callable_integrand_functions(self):
    
        num_of_dims = self.phys_space.num_of_dims
        num_of_nds = self.template_el.num_of_nds

        # Shape functions, used for both trial & test functions
        shape_fns_exprs = aux.flatten_nested_list(self.template_el.ref_shape_fns_exprs, num_of_dims, 0)
        # Other necessary symbols for lambafication
        vec_crd_syms = list(self.template_el.ref_space.vec_crd_syms)
        nds_phys_vec_crds_syms = aux.flatten_nested_list(self.template_el.nds_phys_vec_crds_syms)
        bdry_norm_vec_syms = list(self.placeholder_bdry_norm_vec_expr)
        
        # Storage for callable integrand functions
        vol_integrand_funcs = []
        bdry_integrand_funcs = []
        src_integrand_funcs = []
        
        # Template weak-form PDE components
        weak_op_vol_part = 
        weak_op_bdry_part = 
        weak_src_expr = 

        for curr_nd_i in range(num_of_nds):
            curr_test_fn_expr = shape_fns_exprs[curr_nd_i]
            curr_trial_fn_expr = curr_test_fn_expr # Galerkin method
            
            # Create integrand expressions
            curr_vol_integrand_expr = 
            curr_bdry_integrand_expr = 
            curr_src_integrand_expr = 
            
            # Create callable functions
            vol_integrand_funcs.append(
                sp.lambdify(
                    nds_phys_vec_crds_syms + vec_crd_syms + bdry_norm_vec_syms,
                    curr_vol_integrand_expr,
                    'numpy'
                )
            )
            
            bdry_integrand_funcs.append(
                sp.lambdify(
                    nds_phys_vec_crds_syms + vec_crd_syms + bdry_norm_vec_syms, 
                    curr_bdry_integrand_expr,
                    'numpy'
                )
            )
            
            src_integrand_funcs.append(
                sp.lambdify(
                    nds_phys_vec_crds_syms + vec_crd_syms,
                    curr_src_integrand_expr, 
                    'numpy'
                )
            )
        
        return vol_integrand_funcs, bdry_integrand_funcs, src_integrand_funcs