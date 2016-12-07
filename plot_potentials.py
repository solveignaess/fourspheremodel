import numpy as np
import matplotlib.pyplot as plt

from parameters import *
# Homogeneous infinite medium
I = 1.
d = .1
phi_0 = I*d / (4*np.pi*sigma_brain*(scalp_rad**2))

n = np.arange(0, 200)
coef = n*((dipole_loc / scalp_rad)**(n+1))
Lprod = np.polynomial.legendre.Legendre(coef)
phi_sphere = I*d*Lprod(np.cos(np.deg2rad(theta))) / (4*np.pi*sigma_brain*(dipole_loc**2))

#phi_20_simple = (34.1*np.exp(-0.15*theta) +1.29 - (0.123*theta) + 0.00164*(theta**2))*phi_0
#phi_40 = (27.4*np.exp(-0.10*theta) -5.49 + (0.203*theta) - 0.00234*(theta**2))*phi_0
#phi_80 = (13.4*np.exp(-0.10*theta) - 0.155 - (0.0135*theta))*phi_0

phi_20 = np.load('phi_20.npy')
phi_20_c = np.load('phi_20_c.npy')

phi_40 = np.load('phi_40.npy')
phi_40_c = np.load('phi_40_c.npy')

phi_80 = np.load('phi_80.npy')
phi_80_c = np.load('phi_80_c.npy')

num_20 = np.load('4Shell_FEM_20.npy')
num_40 = np.load('4Shell_FEM_40.npy')
num_80 = np.load('4Shell_FEM_80.npy')

plt.subplot(111)
plt.plot(theta, phi_20, 'k', label='Nunez 20')
plt.plot(theta, phi_20_c, 'b', label='Correction')
plt.plot(theta, num_20,  'r.', label='FEM 20')
plt.xlabel('Polar angle (degrees)')
plt.ylabel('Potential $\mu V$')
plt.title('$\sigma_{skull} =sigma_{brain} / 20$')
plt.legend()
plt.savefig('Potentials_20.png')

plt.clf()
plt.subplot(111)
plt.plot(theta, phi_40, 'k', label='Nunez 40')
plt.plot(theta, phi_40_c, 'b', label='Correction')
plt.plot(theta, num_40,  'r.', label='FEM 40')
plt.xlabel('Polar angle (degrees)')
plt.ylabel('Potential $\mu V$')
plt.title('$\sigma_{skull} =sigma_{brain} / 40$')
plt.legend()
plt.savefig('Potentials_40.png')

plt.clf()
plt.subplot(111)
plt.plot(theta, phi_80, 'k', label='Nunez 80')
plt.plot(theta, phi_80_c, 'b', label='Correction')
plt.plot(theta, num_80,  'r.', label='FEM 80')
plt.xlabel('Polar angle (degrees)')
plt.ylabel('Potential $\mu V$')
plt.title('$\sigma_{skull} =sigma_{brain} / 80$')
plt.legend()
plt.savefig('Potentials_80.png')
#plt.show()
