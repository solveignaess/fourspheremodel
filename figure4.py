# -*- coding: utf-8 -*-
import os
import numpy as np
import matplotlib.pyplot as plt
from plotting_convention import mark_subplots, simplify_axes
import parameters as params
import sys
reload(sys)
sys.setdefaultencoding('UTF-8')

# Homogeneous sphere, Nunez 2006, Eq. 6.7
I = 1.
d = .1
theta = params.theta.reshape(180, 180)[:, 0][0:90]
k = 1.0 #0.2  # multiply I by 200000 and convert from nV to mV

f = params.dipole_loc / params.scalp_rad
K = 1.+ (f**2) - (2*f*np.cos(params.theta_r))
frst_trm = (2.*(np.cos(params.theta_r) - f)) / (K**1.5)
scnd_trm = (1./f)*((1/(K)**0.5) - 1.)
prod = frst_trm + scnd_trm
phi_sphere = I*d*prod / (4*np.pi*params.sigma_brain*(params.scalp_rad**2))*k
phi_sphere = phi_sphere.reshape(180, 180)[:, 0][0:90]

nunsri06 = np.load(os.path.join('results', 'Analytical_NunSri06_rad.npz'))
sri98 = np.load(os.path.join('results', 'Analytical_Sri98_rad.npz'))
analytical = np.load(os.path.join('results', 'Analytical_rad.npz'))

phi_nunsri06 = nunsri06['phi_lim'].reshape(180, 180)[:, 0][0:90] * k
phi_sri98 = sri98['phi_lim'].reshape(180, 180)[:, 0][0:90] * k
phi_correct = analytical['phi_lim'].reshape(180, 180)[:, 0][0:90] * k

fig = plt.figure(figsize=[5, 4])
fig.subplots_adjust(left=0.2, bottom=0.14)
plt.subplot(111)
plt.xlim([0, 50])
plt.ylim([-0.2, 5.])

# plt.ylim([-0.5, 1.+1e-10]) # , xlim=[0, 50], ylim=[-0.2, 5])
# plt.plot(theta, phi_20, 'k', label='Nunez 20')

plt.plot(theta, phi_sphere, 'r', label='Homogeneous sphere (Eq. 9)')
plt.plot(theta, phi_nunsri06, 'g+', label='Nunez & Srinivasan (2006)')
plt.plot(theta, phi_sri98, 'k*', label='Srinivasan (1998)')
plt.plot(theta, phi_correct, 'b.', label='Present results - analytical')

plt.xlabel('Polar angle (degrees)', fontsize=14)
plt.ylabel('Potential (mV)', fontsize=14.)
plt.tick_params(labelsize=15.)
plt.title('$\sigma_{\mathrm{skull}} = \sigma_{\mathrm{brain}} = \sigma_{\mathrm{csf}} = \sigma_{\mathrm{scalp}}$', fontsize=19)
plt.legend(frameon=False, bbox_to_anchor=(1, 0.9), fontsize=11)
simplify_axes(fig.axes)

plt.savefig(os.path.join('results', 'sanity_200000.png'), dpi=150)
#plt.show()
