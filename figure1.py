import os
import numpy as np
import matplotlib.pyplot as plt

from parameters_wm import *

# Homogeneous infinite medium
I = 1.
d = .1
phi_0 = I*d / (4*np.pi*sigma_brain*(scalp_rad**2))

n = np.arange(0, 200)
coef = n*((dipole_loc / scalp_rad)**(n+1))
Lprod = np.polynomial.legendre.Legendre(coef)
phi_sphere = I*d*Lprod(np.cos(np.deg2rad(theta))) / (4*np.pi*sigma_brain*(dipole_loc**2))

phi_20_fit = (34.1*np.exp(-0.15*theta) +1.29 - (0.123*theta) + 0.00164*(theta**2))*phi_0
phi_40_fit = (27.4*np.exp(-0.10*theta) -5.49 + (0.203*theta) - 0.00234*(theta**2))*phi_0
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

plt.figure(figsize=(15,6))

plt.subplot(131)
plt.plot(theta, phi_20, 'k', label='Srinivasan-06')
plt.plot(theta, phi_20_98, 'g', label='Srinivasan-98')
plt.plot(theta[:50], phi_20_fit[:50], 'm+', label='Srinivasan-06-Fit')
plt.plot(theta, phi_20_c, 'b', label='Correction')
plt.plot(theta, num_20,  'r.', label='FEM')
plt.xlabel('Polar angle (degrees)')
plt.ylabel('Potential $\mu V$')
plt.title('$\sigma_{skull} = \sigma_{brain} / 20$')
ax = plt.gca()
ax.set_ylim([-0.2, 1.2])
#plt.legend()
#plt.savefig(os.path.join('results', 'Potentials_20.png'))

#plt.clf()
plt.subplot(132)
plt.plot(theta, phi_40, 'k', label='Srinivasan-06')
plt.plot(theta, phi_40_98, 'g', label='Srinivasan-98')
plt.plot(theta[:50], phi_40_fit[:50], 'm+', label='Srinivasan-06-Fit')
plt.plot(theta, phi_40_c, 'b', label='Correction')
plt.plot(theta, num_40,  'r.', label='FEM')
plt.xlabel('Polar angle (degrees)')
#plt.ylabel('Potential $\mu V$')
plt.title('$\sigma_{skull} = \sigma_{brain} / 40$')
ax = plt.gca()
ax.set_ylim([-0.2, 1.2])
#plt.legend()
#plt.savefig(os.path.join('results', 'Potentials_40.png'))

#plt.clf()
plt.subplot(133)
plt.plot(theta, phi_80, 'k', label='Srinivasan-06')
plt.plot(theta, phi_80_98, 'g', label='Srinivasan-98')
plt.plot(theta[:50], phi_80_fit[:50], 'm+', label='Srinivasan-06-Fit')
plt.plot(theta, phi_80_c, 'b', label='Correction')
plt.plot(theta, num_80,  'r.', label='FEM')
plt.xlabel('Polar angle (degrees)')
#plt.ylabel('Potential $\mu V$')
plt.title('$\sigma_{skull} = \sigma_{brain} / 80$')
ax = plt.gca()
ax.set_ylim([-0.2, 1.2])
plt.legend()
plt.savefig(os.path.join('results', 'Potentials.png'), dpi=300)
#plt.show()
