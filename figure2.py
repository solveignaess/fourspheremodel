from __future__ import division
import numpy as np
from CalcPotential4Sphere import CalcPotential4Sphere
import matplotlib.pyplot as plt
from matplotlib import colors
from mpl_toolkits.mplot3d import Axes3D

'''Plot scalp potentials calculated analytically with the four-sphere model,
and numerically with FEM.'''

plt.close('all')

# load scalp potentials from eeg simulations
eeg_rad = np.load('./results/eeg_radial.npz')
eeg_tan = np.load('./results/eeg_tangential.npz')
eeg_mix = np.load('./results/eeg_mix.npz')
files = [eeg_rad, eeg_tan , eeg_mix]

# create array with 2D electrode positions for plotting geometry
radii = [79000., 80000., 85000., 90000.]
r = 90000.
xs = np.linspace(-r, r, 1001)
zs = np.linspace(-r, r, 1001)
elec_pos = []
for zpos in zs:
    for xpos in xs:
        elec_pos.append([xpos, 0., zpos])
elec_pos = np.array(elec_pos)

# mark subdomains for plotting geometry
subdomain_markers = []
for pos in elec_pos:
    rad = np.sqrt(pos[0]**2 + pos[1]**2 + pos[2]**2)
    # print r
    if rad < radii[0]:  #  or (r-79000.) < 1e-16:
        subdomain_markers.append(0.1)
    elif rad < radii[1]: #  or (r-80000.) < 1e-16:
        subdomain_markers.append(0.2)
    elif rad < radii[2]: #  or (r-90000.) < 1e-16:
        subdomain_markers.append(0.3)
    elif rad < radii[3]:  # or (r-100000.) < 1e-16:
        subdomain_markers.append(0.4)
    else:
        subdomain_markers.append(np.nan)
subd_markers = np.array(subdomain_markers).reshape(len(xs), len(zs))

# create array with 3d scalp electrode positions
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
X = (r-1e2) * np.outer(np.cos(u), np.sin(v))
Y = (r-1e2) * np.outer(np.sin(u), np.sin(v))
Z = (r-1e2) * np.outer(np.ones(np.size(u)), np.cos(v))

# compute eeg analytically,
# compute absolute error between analytical and fem calculation
phi_4s_list = [np.zeros(X.shape), np.zeros(X.shape), np.zeros(X.shape)]
error_list = [np.zeros(X.shape), np.zeros(X.shape), np.zeros(X.shape)]

P1s = []
rz1s = []
phi_fem_list = []
for idx, f in enumerate(files):
    param_dict = f['params'].item()
    fem_phi = f['phi']
    phi_fem_list.append(fem_phi)
    charge_pos = param_dict['charge_pos']
    charges = param_dict['charge']
    rz1 = (charge_pos[0] + charge_pos[1])/2
    P1 = np.array([charge_pos[0]*charges[0] + charge_pos[1]*charges[1]])
    P1s.append(P1)
    r = param_dict['radii'][3]
    radii = param_dict['radii']
    sigmas = param_dict['sigmas']
    phi_4s = phi_4s_list[idx]
    error = error_list[idx]
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            sphere_mod = CalcPotential4Sphere(radii, sigmas, np.array([X[i,j], Y[i,j], Z[i,j]]), rz1)
            pot = sphere_mod.calc_potential(P1)
            phi_4s[i,j] = np.nan_to_num(pot)
            error[i,j] = np.abs(phi_4s[i,j] - fem_phi[i,j])
    print 'error max:', np.max(np.abs(error))

# create color lists for plotting of potentials and error
vmax = .1
vmin = -vmax
clr = lambda phi: plt.cm.PRGn((phi - vmin) / (vmax - vmin))
vmax_error = 0.006
vmin_error = 0.
clr_error = lambda error: plt.cm.Blues((error - vmin_error) / (vmax_error - vmin_error))

colors_4s_list = [clr(phi_4s_list[idx]) for idx in range(3)]
colors_fem_list = [clr(phi_fem_list[idx]) for idx in range(3)]
colors_error_list = [clr_error(error_list[idx]) for idx in range(3)]

### plotting ###
plt.close('all')
fig = plt.figure()

ax1 = plt.subplot2grid((3,4),(0,0))
ax2 = plt.subplot2grid((3,4),(0,1), projection='3d')
ax3 = plt.subplot2grid((3,4),(0,2), projection='3d')
ax4 = plt.subplot2grid((3,4),(0,3), projection='3d')
ax5 = plt.subplot2grid((3,4),(1,0))
ax6 = plt.subplot2grid((3,4),(1,1), projection='3d')
ax7 = plt.subplot2grid((3,4),(1,2), projection='3d')
ax8 = plt.subplot2grid((3,4),(1,3), projection='3d')
ax9 = plt.subplot2grid((3,4),(2,0))
ax10 = plt.subplot2grid((3,4),(2,1), projection='3d')
ax11 = plt.subplot2grid((3,4),(2,2), projection='3d')
ax12 = plt.subplot2grid((3,4),(2,3), projection='3d')
ax2.set_title('4-sphere', y = 1.1)
ax3.set_title('FEM', y = 1.1)
ax4.set_title('Error', y = 1.1)

# plot geometry
for ax, P1 in zip([ax1, ax5, ax9], P1s):
    im = ax.imshow(subd_markers.reshape(len(xs), len(zs)).T,
                    origin='lower', interpolation='nearest', extent = [-90000., 90000., -90000., 90000.])
    cmap = colors.ListedColormap([ '#999999', '#83caff', '#e6e6e6', '#ffcc99'])
    im.set_cmap(cmap)
    ax.add_patch(plt.Circle((0, 0), radius = radii[0], color = 'gray', fill=False, lw = .3))
    ax.add_patch(plt.Circle((0, 0), radius = radii[1], color = 'gray', fill=False, lw = .3))
    ax.add_patch(plt.Circle((0, 0), radius = radii[2], color = 'gray', fill=False, lw = .3))
    ax.add_patch(plt.Circle((0,0), radius = radii[3], color = 'gray', fill=False, lw = .5))
    ax.arrow(-100000, -100000, 0, 18000, head_width=4000, head_length=4000, fc='k', ec='k')
    ax.arrow(-100000, -100000, 18000,0, head_width=4000, head_length=4000, fc='k', ec='k')
    ax.arrow(-100000, -100000, -7800, -7800, head_width=4000, head_length=3300, fc='k', ec='k')
    ax.text(-102800, -70000, 'z', size = 6.)
    ax.text(-70000, -103800, 'x', size = 6.)
    ax.text(-115000, -120000, 'y', size = 6.)
    ax.set_xlim(-120000., 120000.)
    ax.set_ylim(-120000., 120000.)

    arrow = np.sum(P1, axis = 0)
    if ax == ax1:
        px = charge_pos[0][0]
        pz = charge_pos[0][2]-20*1000
    elif ax == ax5:
        px = charge_pos[0][0] + 3000
        pz = charge_pos[0][2] - 10000  # + 2250
    else:
        px = charge_pos[0][0] + 1590
        pz = charge_pos[0][2]- 20*1000+1590
    ax.arrow(px, pz,
              11*arrow[0], 11*arrow[2],
              fc='k',
              ec='k',
              width = 170,
              length_includes_head=False)
    ax.axis('off')
# plot analytically calculated potentials
for ax, clrs in zip([ax2, ax6, ax10], colors_4s_list):
    surf1 = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=clrs,
                             linewidth=0, antialiased=False)
# plot FEM-computed potentials
for ax, clrs in zip([ax3, ax7, ax11], colors_fem_list):
    surf2 = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=clrs,
                             linewidth=0, antialiased=False)
# plot errors
for ax, clrs in zip([ax4, ax8, ax12], colors_error_list):
    surf3 = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors = clrs,#cmap=plt.cm.Blues,
                   linewidth=0, antialiased=False)  #, norm=colors.LogNorm(vmin=1e-7, vmax=1e-2))

for ax in [ax2, ax3, ax4, ax6, ax7, ax8, ax10, ax11, ax12]:
    ax.set_aspect('equal')
    ax.axis('off')
    # ax.auto_scale_xyz([-60000., 60000.], [-60000., 60000.], [-60000., 60000.])
    ax.set_xlim3d(-65000, 65000)
    ax.set_ylim3d(-65000, 65000)
    ax.set_zlim3d(3000-65000, 3000+65000)
    ax.view_init(10, 270)

# colorbars
cax1 = fig.add_axes([0.31, 0.07, 0.38, 0.01])
m = plt.cm.ScalarMappable(cmap=plt.cm.PRGn)
m.set_array(phi_4s_list[1])
cbar1 = fig.colorbar(m, cax=cax1, format='%3.3f', extend = 'both', orientation='horizontal')
cbar1.outline.set_visible(False)
ticks = np.linspace(-0.4, 0.4, 9)
plt.xticks(ticks, [str(tick) for tick in ticks], rotation = 40)
cbar1.set_label('Potential (nV)', labelpad=5.2)

cax2 = fig.add_axes([0.79, 0.07, 0.15, 0.01])
m = plt.cm.ScalarMappable(cmap=plt.cm.Blues)
ticks2 = np.linspace(0.001, 0.005, 9)
m.set_array(error_list[0])
cbar2 = fig.colorbar(m, cax=cax2, format='%3.6f', extend='max', orientation='horizontal')
cbar2.outline.set_visible(False)
cax2.set_xticks(ticks2)
cax2.set_xticklabels(['0.001', '', '0.002', '', '0.003', '', '0.004', '', '0.005'], rotation = 40)
cbar2.set_label('Potential (nV)', labelpad=.1)

fig.subplots_adjust(wspace = 0., hspace=-0.01, right = 0.98, left = 0.02)

# label figures
fig.text(0.03, .87, 'A',
        horizontalalignment='center',
        verticalalignment='center',
        fontweight='demibold',
        fontsize=12)
fig.text(0.27, .87, 'B',
        horizontalalignment='center',
        verticalalignment='center',
        fontweight='demibold',
        fontsize=12)
fig.text(0.51, .87, 'C',
        horizontalalignment='center',
        verticalalignment='center',
        fontweight='demibold',
        fontsize=12)
fig.text(0.75, .87, 'D',
        horizontalalignment='center',
        verticalalignment='center',
        fontweight='demibold',
        fontsize=12)
fig.text(0.03, .61, 'E',
        horizontalalignment='center',
        verticalalignment='center',
        fontweight='demibold',
        fontsize=12)
fig.text(0.27, .61, 'F',
        horizontalalignment='center',
        verticalalignment='center',
        fontweight='demibold',
        fontsize=12)
fig.text(0.51, .61, 'G',
        horizontalalignment='center',
        verticalalignment='center',
        fontweight='demibold',
        fontsize=12)
fig.text(0.75, .61, 'H',
        horizontalalignment='center',
        verticalalignment='center',
        fontweight='demibold',
        fontsize=12)
fig.text(0.03, .34, 'I',
        horizontalalignment='center',
        verticalalignment='center',
        fontweight='demibold',
        fontsize=12)
fig.text(0.27, .34, 'J',
        horizontalalignment='center',
        verticalalignment='center',
        fontweight='demibold',
        fontsize=12)
fig.text(0.51, .34, 'K',
        horizontalalignment='center',
        verticalalignment='center',
        fontweight='demibold',
        fontsize=12)
fig.text(0.75, .34, 'L',
        horizontalalignment='center',
        verticalalignment='center',
        fontweight='demibold',
        fontsize=12)


fig.set_size_inches(9., 6.)
plt.savefig('./results/eeg_fig2.pdf', dpi=600., bbox_inches='tight')
