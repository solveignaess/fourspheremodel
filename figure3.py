import os
import numpy as np
import matplotlib.pyplot as plt
from plotting_convention import mark_subplots, simplify_axes
from parameters_wm import *

# Homogeneous infinite medium
I = 200000.
d = .1
phi_0 = I*d / (4*np.pi*sigma_brain*(scalp_rad**2))

k = 1e-6  # convert from nV to mV

n = np.arange(0, 200)
coef = n*((dipole_loc / scalp_rad)**(n+1))
Lprod = np.polynomial.legendre.Legendre(coef)
phi_sphere = I*d*Lprod(np.cos(np.deg2rad(theta))) / (4*np.pi*sigma_brain*(dipole_loc**2))

phi_20_fit = (34.1*np.exp(-0.15*theta) + 1.29 - (0.123*theta) + 0.00164*(theta**2))*phi_0*k
phi_40_fit = (27.4*np.exp(-0.10*theta) - 5.49 + (0.203*theta) - 0.00234*(theta**2))*phi_0*k
phi_80_fit = (13.4*np.exp(-0.10*theta) - 0.155 - (0.0135*theta))*phi_0*k

phi_20 = np.load(os.path.join('results', 'phi_20_200000.npy'))*k
phi_20_98 = np.load(os.path.join('results', 'phi_20_98_200000.npy'))*k
phi_20_c = np.load(os.path.join('results', 'phi_20_c_200000.npy'))*k

phi_40 = np.load(os.path.join('results', 'phi_40_200000.npy'))*k
phi_40_98 = np.load(os.path.join('results', 'phi_40_98_200000.npy'))*k
phi_40_c = np.load(os.path.join('results', 'phi_40_c_200000.npy'))*k

phi_80 = np.load(os.path.join('results', 'phi_80_200000.npy'))*k
phi_80_98 = np.load(os.path.join('results', 'phi_80_98_200000.npy'))*k
phi_80_c = np.load(os.path.join('results', 'phi_80_c_200000.npy'))*k

num_20 = np.load(os.path.join('results', '4Shell_FEM_20_wm_high_res2.npy'))*k
num_40 = np.load(os.path.join('results', '4Shell_FEM_40_wm_high_res2.npy'))*k
num_80 = np.load(os.path.join('results', '4Shell_FEM_80_wm_high_res2.npy'))*k

fig = plt.figure(figsize=(12, 4))
fig.subplots_adjust(wspace = 0.28, bottom=0.15, left=0.08, right=0.99)

ylim = [-0.050000000000000003, 0.25000000000000006]
xlim = [0., 80.00000000000000003]
ax = plt.subplot(131, ylim=ylim, xlim=xlim)
ax.set_title(r'$\sigma_{\mathrm{skull}} = \sigma_{\mathrm{brain}} / 20$', fontsize=17)
ax.set_ylabel('Potential (mV)', fontsize=14)
plt.plot(theta, phi_20, 'k', lw=2)
plt.plot(theta, phi_20_98, 'g', lw=2)
plt.plot(theta[:50], phi_20_fit[:50], 'm+',  lw=2)
plt.plot(theta, phi_20_c, 'b',  lw=2)
plt.plot(theta, num_20,  'r.', lw=2)
ax.set_xticks(ax.get_xticks()[::2])
ax.tick_params(labelsize=15.)



ax = plt.subplot(132, ylim=ylim, xlim=xlim)
ax.set_xlabel('Polar angle (degrees)', fontsize=14)
ax.set_title(r'$\sigma_{\mathrm{skull}} = \sigma_{\mathrm{brain}} / 40$', fontsize=17)
plt.plot(theta, phi_40, 'k', lw=2)
plt.plot(theta, phi_40_98, 'g', lw=2)
plt.plot(theta[:50], phi_40_fit[:50], 'm+', lw=2)
plt.plot(theta, phi_40_c, 'b', lw=2)
plt.plot(theta, num_40,  'r.', lw=2)
ax.set_xticks(ax.get_xticks()[::2])
ax.tick_params(labelsize=15.)

ax = plt.subplot(133, ylim=ylim, xlim=xlim)
ax.set_title(r'$\sigma_{\mathrm{skull}} = \sigma_{\mathrm{brain}} / 80$', fontsize=17)
l1, = plt.plot(theta, phi_80_98, 'g', lw=2)
l2, = plt.plot(theta, phi_80, 'k', lw=2)
l3, = plt.plot(theta[:50], phi_80_fit[:50], 'm+', lw=2)
l4, = plt.plot(theta, phi_80_c, 'b', lw=2)
l5, = plt.plot(theta, num_80,  'r.', lw=2)
ax.set_xticks(ax.get_xticks()[::2])
ax.tick_params(labelsize=15.)

lines = [l1, l2, l3, l4, l5]
line_names = ['Srinivasan (1998)', 'Nunez & Srinivasan (2006)',
              'Fit - Nunez & Srinivasan (2006) ', 'Present results - analytical', 'Present results - FEM']
fig.legend(lines, line_names, frameon=False, ncol=1, bbox_to_anchor=[1.01, 0.85], fontsize=12)

simplify_axes(fig.axes)
mark_subplots(fig.axes, ypos=1.07, xpos=-0.)
plt.savefig(os.path.join('results', 'Potentials3.png'), dpi=150)
#plt.show()
