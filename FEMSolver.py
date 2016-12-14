"""
Solve PDE: -div(sigma div u) = f, with du/dn = g on boundary.

"""
from __future__ import division

from dolfin import *
from mshr import *
import numpy as np


class FEMSolver:
    def __init__(self, param_dict, mesh=None, subdomains=None):
        self.mesh = mesh
        self.subdomains = subdomains
        self.sigmas = [Constant(s) for s in param_dict['sigmas']]
        self.charge_pos = param_dict['charge_pos']
        self.charge = param_dict['charge']

    def fem_pots(self, poly_order=1):
        """Calculate potentials in 4-sphere head model

        Parameters
        __________
        mesh : xml-file - head mesh.
              Preferably 4 concentric spheres with smooth boundaries
        charge_pos : list/ ndarray [x, y, z] in [micro m]
                     Current source locations in brain.
        charge : list/ ndarray [nA]
                 Current strengths, same order as charge_pos

        Returns
        _______
        phi : array []
        """
        parameters["krylov_solver"]["relative_tolerance"] = 1e-10
        V = FunctionSpace(self.mesh, "CG", poly_order)
        u = TrialFunction(V)
        v = TestFunction(V)
        dx = Measure("dx")[self.subdomains]

        whitemattervol = 32
        graymattervol = 64
        csfvol = 96
        skullvol = 128
        scalpvol = 160

        a = (inner(self.sigmas[0] * grad(u), grad(v)) * dx(whitemattervol) +
            inner(self.sigmas[0] * grad(u), grad(v)) * dx(graymattervol) +
            inner(self.sigmas[1] * grad(u), grad(v)) * dx(csfvol) +
            inner(self.sigmas[2] * grad(u), grad(v)) * dx(skullvol) +
            inner(self.sigmas[3] * grad(u), grad(v)) * dx(scalpvol))

        L = Constant(0)*v*dx

        A, b = assemble_system(a, L)

        # set point sources as boundary conditions
        for idx, pos in enumerate(self.charge_pos):
            point = Point(pos[0], pos[1], pos[2])
            delta = PointSource(V, point, self.charge[idx])
            delta.apply(b)

        phi = Function(V)

        solver = KrylovSolver("cg", "ilu")
        solver.parameters["maximum_iterations"] = 1000
        solver.parameters["absolute_tolerance"] = 1E-8
        solver.parameters["monitor_convergence"] = True

        info(solver.parameters, True)
        set_log_level(PROGRESS)
        solver.solve(A, phi.vector(), b)

        k1 = 1E2 # from 10**-7V to nV
        phi = phi*k1
        return phi
