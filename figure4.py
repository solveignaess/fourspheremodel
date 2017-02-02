# -*- coding: utf-8 -*-
import os
import numpy as np
import matplotlib.pyplot as plt
from plotting_convention import mark_subplots, simplify_axes
from parameters_wm import *
import sys
reload(sys)
sys.setdefaultencoding('UTF-8')

# Homogeneous sphere, Nunez 2006, Eq. 6.7
I = 1.
d = .1

k = 0.2  # multiply I by 200000 and convert from nV to mV

f = dipole_loc / scalp_rad
K = 1.+ (f**2) - (2*f*np.cos(np.deg2rad(theta)))
frst_trm = (2.*(np.cos(np.deg2rad(theta)) - f)) / (K**1.5)
scnd_trm = (1./f)*((1/(K)**0.5) - 1.)
prod = frst_trm + scnd_trm
phi_sphere = I*d*prod / (4*np.pi*sigma_brain*(scalp_rad**2))*k

#phi_20_simple = (34.1*np.exp(-0.15*theta) +1.29 - (0.123*theta) + 0.00164*(theta**2))*phi_0
#phi_40 = (27.4*np.exp(-0.10*theta) -5.49 + (0.203*theta) - 0.00234*(theta**2))*phi_0
#phi_80 = (13.4*np.exp(-0.10*theta) - 0.155 - (0.0135*theta))*phi_0

phi_20_sanity = np.load(os.path.join('results', 'phi_sanity.npy'))*k
phi_20_sanity_98 = np.load(os.path.join('results', 'phi_sanity_98.npy'))*k
phi_20_c_sanity = np.load(os.path.join('results', 'phi_c_sanity.npy'))*k

fig = plt.figure(figsize=[5, 4])
fig.subplots_adjust(left=0.2, bottom=0.14)
plt.subplot(111)
plt.ylim([-0.05, 1.+1e-10]) # , xlim=[0, 50], ylim=[-0.2, 5])
#plt.plot(theta, phi_20, 'k', label='Nunez 20')

plt.plot(theta, phi_sphere,  'r', label='Homogeneous sphere (Eq. 9)')
plt.plot(theta, phi_20_sanity, 'g+', label='Nunez & Srinivasan (2006)')
plt.plot(theta, phi_20_sanity_98, 'k*', label='Srinivasan (1998)')
plt.plot(theta, phi_20_c_sanity, 'b.', label='Present results - analytical')

plt.xlabel('Polar angle (degrees)', fontsize=14)
plt.ylabel('Potential (mV)', fontsize=14.)
plt.tick_params(labelsize=15.)
plt.title('$\sigma_{\mathrm{skull}} = \sigma_{\mathrm{brain}} = \sigma_{\mathrm{csf}} = \sigma_{\mathrm{scalp}}$', fontsize=19)
plt.legend(frameon=False, bbox_to_anchor=(1, 0.9), fontsize=11)
simplify_axes(fig.axes)
plt.savefig(os.path.join('results', 'sanity_200000.png'), dpi=300)

#plt.show()
