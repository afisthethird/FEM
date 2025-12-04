# Libraries
import numpy as np
# Scripts
import PYTHON.space as spce
import PYTHON.equation as eqn
import PYTHON.mesh as msh


class Simulation():
    def __init__(
        self,
        name                : str,
        fem_eqn             : eqn.FiniteElementMethodEquation,
        phys_dom            : spce.Domain,
        nums_of_els_per_dim : list[int]
        ):

        self.name = name
        self.fem_eqn = fem_eqn
        # Mesh
        self.mesh = msh.Mesh(phys_dom, self.fem_eqn.template_el, nums_of_els_per_dim)