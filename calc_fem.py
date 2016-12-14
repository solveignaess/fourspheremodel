from __future__ import division
from dolfin import *
import numpy as np
from FEMSolver import FEMSolver


r1 = 7.9; r2 = 8.0; r3 = 8.5; r4 = 9.0  # [cm]
res1 = 170; res2 = 180; res3 = 190; res4 = 200  # res1 = 10; res2 = 250; res3 = 150; res4 = 150
charge_pos = np.array([[0., 0., 7.79], [0., 0., 7.81]], dtype=float)  # [cm]
charge = [10, -10]  # [nA]

param_dict = {'radii': [r1, r2, r3, r4],
              'res_spheres': [res1, res2, res3, res4],
              'eps': 1e-5,
              'charge_pos': charge_pos,
              'charge': charge}
param_dict['sigmas'] = [0.3, 1.5, 0.015, 0.3]

mesh = Mesh('./mesh/sphere_4_wm.xml')
subdomains = MeshFunction("size_t", mesh, './mesh/sphere_4_wm_physical_region.xml')
fem_solver = FEMSolver(param_dict, mesh, subdomains)
phi = fem_solver.fem_pots(poly_order = 2)

# only surface values
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)

X = (r4-1e-2) * np.outer(np.cos(u), np.sin(v))
Y = (r4-1e-2) * np.outer(np.sin(u), np.sin(v))
Z = (r4-1e-2) * np.outer(np.ones(np.size(u)), np.cos(v))

values = np.zeros(X.shape)
for i in range(X.shape[0]):
    for j in range(X.shape[1]):
        val = phi(np.array([X[i,j], Y[i,j], Z[i,j]]))
        values[i,j] = val

# update dictionary radius values from cm to mu m
param_dict['radii'] = [s*1e4 for s in param_dict['radii']]
param_dict['charge_pos'] = [p*1e4 for p in param_dict['charge_pos']]
# now the stored values are in units mu m, nA and nV
np.savez('./results/chaitanya_wm_sigma20_poly2_eeg.npz', params = param_dict, phi = values)
