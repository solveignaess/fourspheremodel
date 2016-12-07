import numpy as np

# All numbers in cm
dipole_loc = 7.8
brain_rad = 7.9
csftop_rad = 8.
skull_rad = 8.5
scalp_rad = 9.

sigma_brain = 1. / 300.  # S / cm
sigma_scalp = sigma_brain
sigma_csf = 5 * sigma_brain
sigma_skull20 = sigma_brain / 20.
sigma_skull40 = sigma_brain / 40.
sigma_skull80 = sigma_brain / 80.

# from gmsh
brainvol = 32
csfvol = 64
skullvol = 96
scalpvol = 128
airvol = 160
airsurf = 158

# measument points
theta = np.arange(1, 90)
phi_angle = 0
rad_tol = 1e-2
x_points = (scalp_rad - rad_tol) * np.sin(np.deg2rad(theta)) * np.cos(np.deg2rad(phi_angle))
y_points = (scalp_rad - rad_tol) * np.sin(np.deg2rad(theta)) * np.sin(np.deg2rad(phi_angle))
z_points = (scalp_rad - rad_tol) * np.cos(np.deg2rad(theta))

ele_coords = np.vstack((x_points, y_points, z_points)).T

# dipole location
src_pos = [0., 0., 7.85]
snk_pos = [0., 0., 7.75]
# src_pos = [0., 0., 0.05]
# snk_pos = [0., 0., -0.05]
