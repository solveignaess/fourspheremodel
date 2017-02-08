import os
import numpy as np
from dolfin import *
import parameters as params

def extract_pots(phi, positions):
    compt_values = np.zeros(positions.shape[0])
    for ii in range(positions.shape[0]):
        compt_values[ii] = phi(positions[ii, :])
    return compt_values

def main_4shell_fem(mesh, subdomains, boundaries, skull_cond, src_pos, snk_pos):
    sigma_B = Constant(params.sigma_brain)
    sigma_Sc = Constant(params.sigma_scalp)
    sigma_C = Constant(params.sigma_csf)
    sigma_Sk = Constant(skull_cond)
    sigma_A = Constant(0.)

    V = FunctionSpace(mesh, "CG", 2)
    v = TestFunction(V)
    u = TrialFunction(V)

    phi = Function(V)
    dx = Measure("dx")[subdomains]
    ds = Measure("ds")[boundaries]
    a = inner(sigma_B * grad(u), grad(v))*dx(params.whitemattervol) + \
        inner(sigma_B * grad(u), grad(v))*dx(params.graymattervol) + \
        inner(sigma_Sc * grad(u), grad(v))*dx(params.scalpvol) + \
        inner(sigma_C * grad(u), grad(v))*dx(params.csfvol) + \
        inner(sigma_Sk * grad(u), grad(v))*dx(params.skullvol)
    L = Constant(0)*v*dx
    A = assemble(a)
    b = assemble(L)

    x_pos, y_pos, z_pos = src_pos
    point = Point(x_pos, y_pos, z_pos)
    delta = PointSource(V, point, 1.)
    delta.apply(b)

    x_pos, y_pos, z_pos = snk_pos
    point1 = Point(x_pos, y_pos, z_pos)
    delta1 = PointSource(V, point1, -1.)
    delta1.apply(b)

    solver = KrylovSolver("cg", "ilu")
    solver.parameters["maximum_iterations"] = 1000
    solver.parameters["absolute_tolerance"] = 1E-8
    solver.parameters["monitor_convergence"] = True

    info(solver.parameters, True)
    set_log_level(PROGRESS)
    solver.solve(A, phi.vector(), b)

    ele_pos_list = params.ele_coords
    vals = extract_pots(phi, ele_pos_list)
    # np.save(os.path.join('results', save_as), vals)
    return vals

if __name__ == '__main__':
    print 'Loading meshes'
    mesh = Mesh(os.path.join("mesh", "sphere_4.xml"))
    subdomains = MeshFunction("size_t", mesh, os.path.join("mesh", "sphere_4_physical_region.xml"))
    boundaries = MeshFunction("size_t", mesh, os.path.join("mesh", "sphere_4_facet_region.xml"))
    for dipole in params.dipole_list:
        print 'Now computing FEM for dipole: ', dipole['name']
        src_pos = dipole['src_pos']
        snk_pos = dipole['snk_pos']
        print 'Done loading meshes'
        fem_20 = main_4shell_fem(mesh, subdomains, boundaries,
                                 params.sigma_skull20, src_pos, snk_pos)
        print 'Done 4Shell-FEM-20'
        fem_40 = main_4shell_fem(mesh, subdomains, boundaries,
                                 params.sigma_skull40, src_pos, snk_pos)
        print 'Done 4Shell-FEM-40'
        fem_80 = main_4shell_fem(mesh, subdomains, boundaries,
                                 params.sigma_skull80, src_pos, snk_pos)
        print 'Done 4Shell-FEM-80'
        f = open(os.path.join('results',
                              'Numerical_' + dipole['name'] + '.npz'), 'w')
        np.savez(f, fem_20=fem_20, fem_40=fem_40, fem_80=fem_80)
        f.close()
