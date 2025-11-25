# Libraries
import sympy as sp
import numpy as np
# Scripts
from Code.types import *
from Code.symbolic.space import PhysicalSpace
from Code.symbolic.math import SymbolicMathematicalObject, Expression, Argument, Operator
from Code.symbolic.math import Derivative, Gradient, Divergence, Laplacian
from Code.symbolic.geometry import Boundary, Domain


class GoverningEquation():
    
    def __init__(
        self,
        name : NameType,
        op : Operator,
        src : Expression,
        host_spce : PhysicalSpace
        ):

        self.name = name
        self.op = op
        self.src = src
        self.host_spce = host_spce
        self.host_spce.create_reference_space()
    
    def construct_weak_form_integrands(
        self,
        wghtng_plchldr : Argument = Argument('w')
        ):
        """
        Constructs the weak-form integrands for the governing equation's LHS (operator on the unknown function) and RHS (source expression) using a placeholder for the weighting functions.
        
        The weighting functions serve to enforce equality at a local level: `∫((L(u)-f)·w)dΩ = 0`.
        """

        # Bookkeeping definitions
        strong_form_op = self.op
        strong_form_src = self.src

        # Conversion
        self.weak_form_op_intgrnd = strong_form_op * wghtng_plchldr
        self.weak_form_src_intgrnd = strong_form_src * wghtng_plchldr
    
    def perform_integration_by_parts(
        self,
        target_deriv_order : NumericIntegerValueType = 1    
        ):
        """
        Performs symbolic (unevaluated) integration by parts on the weak-form operator integrand to reduce derivatives of order > `target_deriv_order` to derivatives of order = `target_deriv_order`.
        
        For a scalar-valued unknown function `u(x): ℝᵈ ⟼ ℝ`, where `x ∈ ℝᵈ`: `∫(∇²u(x)·v(x))dΩ = -∫(∇u·∇v)dΩ + ∫(∇u·vŜ)dS`, where `S = ∂Ω`.
        
        - `-∫(∇u·∇v)dΩ` is known as the "volume" term since it encapsulates the local-to-each-element-subdomain behavior with enforced equality.
        
        -   `∫(∇u·vŜ)dS` is known as the "boundary" term because it is only defined at physical domain boundaries. There is physical meaning associated with integration by parts.
        """
        
        return NotImplementedError
        

# Volume term: -∫(∇u·∇v)dΩ
# Boundary term: ∫(∇u·vŜ)dS
# Gradient: ∇f(x) = J⊥·∇f(ξ)
# Divergence: ∇·v(x) = (1/|J|)*∇·(|J|*J·v(ξ))
