import os
import numpy as np
import matplotlib.pyplot as plt
from plotting_convention import mark_subplots, simplify_axes
from parameters_wm import *

# Homogeneous infinite medium
I = 1.
d = .1
phi_0 = I*d / (4*np.pi*sigma_brain*(scalp_rad**2))

n = np.arange(0, 200)
coef = n*((dipole_loc / scalp_rad)**(n+1))
Lprod = np.polynomial.legendre.Legendre(coef)
phi_sphere = I*d*Lprod(np.cos(np.deg2rad(theta))) / (4*np.pi*sigma_brain*(dipole_loc**2))

phi_20_fit = (34.1*np.exp(-0.15*theta) + 1.29 - (0.123*theta) + 0.00164*(theta**2))*phi_0
phi_40_fit = (27.4*np.exp(-0.10*theta) - 5.49 + (0.203*theta) - 0.00234*(theta**2))*phi_0
phi_80_fit = (13.4*np.exp(-0.10*theta) - 0.155 - (0.0135*theta))*phi_0

phi_20 = np.load(os.path.join('results', 'phi_20.npy'))
phi_20_98 = np.load(os.path.join('results', 'phi_20_98.npy'))
phi_20_c = np.load(os.path.join('results', 'phi_20_c.npy'))

phi_40 = np.load(os.path.join('results', 'phi_40.npy'))
phi_40_98 = np.load(os.path.join('results', 'phi_40_98.npy'))
phi_40_c = np.load(os.path.join('results', 'phi_40_c.npy'))

phi_80 = np.load(os.path.join('results', 'phi_80.npy'))
phi_80_98 = np.load(os.path.join('results', 'phi_80_98.npy'))
phi_80_c = np.load(os.path.join('results', 'phi_80_c.npy'))

num_20 = np.load(os.path.join('results', '4Shell_FEM_20_wm.npy'))
num_40 = np.load(os.path.join('results', '4Shell_FEM_40_wm.npy'))
num_80 = np.load(os.path.join('results', '4Shell_FEM_80_wm.npy'))

fig = plt.figure(figsize=(12, 4))
fig.subplots_adjust(bottom=0.12, left=0.05, right=0.99)

ax = plt.subplot(131, ylim=[-0.1, 1.1], ylabel='Potential ($\mu$V)',
                 title='$\sigma_{skull} = \sigma_{brain} / 20$')
plt.plot(theta, phi_20, 'k', lw=2)
plt.plot(theta, phi_20_98, 'g', lw=2)
plt.plot(theta[:50], phi_20_fit[:50], 'm+',  lw=2)
plt.plot(theta, phi_20_c, 'b',  lw=2)
plt.plot(theta, num_20,  'r.', lw=2)
ax.set_xticks(ax.get_xticks()[::2])


ax = plt.subplot(132, ylim=[-0.1, 1.1], title='$\sigma_{skull} = \sigma_{brain} / 40$',
                 xlabel='Polar angle (degrees)')
plt.plot(theta, phi_40, 'k', lw=2)
plt.plot(theta, phi_40_98, 'g', lw=2)
plt.plot(theta[:50], phi_40_fit[:50], 'm+', lw=2)
plt.plot(theta, phi_40_c, 'b', lw=2)
plt.plot(theta, num_40,  'r.', lw=2)
ax.set_xticks(ax.get_xticks()[::2])

ax = plt.subplot(133, ylim=[-0.1, 1.1], title='$\sigma_{skull} = \sigma_{brain} / 80$')
l1, = plt.plot(theta, phi_80_98, 'g', lw=2)
l2, = plt.plot(theta, phi_80, 'k', lw=2)
l3, = plt.plot(theta[:50], phi_80_fit[:50], 'm+', lw=2)
l4, = plt.plot(theta, phi_80_c, 'b', lw=2)
l5, = plt.plot(theta, num_80,  'r.', lw=2)
ax.set_xticks(ax.get_xticks()[::2])

lines = [l1, l2, l3, l4, l5]
line_names = ['Srinivasan (1998)', 'Nunez & Srinivasan (2006)',
              'Fit - Nunez & Srinivasan (2006) ', 'Correction', 'FEM']
fig.legend(lines, line_names, frameon=False, ncol=1, bbox_to_anchor=[0.99, 0.85], fontsize=11)

simplify_axes(fig.axes)
mark_subplots(fig.axes, ypos=1.05, xpos=-0.)
plt.savefig(os.path.join('results', 'Potentials2.png'), dpi=150)
#plt.show()
