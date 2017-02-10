from __future__ import division
import numpy as np
#from CalcPotential4Sphere import CalcPotential4Sphere
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
import parameters as params

I = 1.0
dipole_list = params.dipole_list
P1s = []
name_list = []
radii = [params.brain_rad, params.csftop_rad,
         params.skull_rad, params.scalp_rad]

for dipole in dipole_list:
    src = np.array(dipole["src_pos"])
    snk = np.array(dipole["snk_pos"])
    P1s.append(np.array((src * I, -1 * snk * I)))
    name_list.append(dipole["name"])

ana_rad = np.load(os.path.join("results", "Analytical_rad.npz"))['phi_20']
ana_tan = np.load(os.path.join("results", "Analytical_tan.npz"))['phi_20']
ana_mix = np.load(os.path.join("results", "Analytical_mix.npz"))['phi_20']

num_rad = np.load(os.path.join("results", "Numerical_rad.npz"))['fem_20']
num_tan = np.load(os.path.join("results", "Numerical_tan.npz"))['fem_20']
num_mix = np.load(os.path.join("results", "Numerical_mix.npz"))['fem_20']

fem_list = [num_rad, num_tan, num_mix]
ana_list = [ana_rad, ana_tan, ana_mix]

fem_list = [ii.reshape(180, 180) for ii in fem_list]
ana_list = [ii.reshape(180, 180) for ii in ana_list]
error_list = [np.abs(ii - jj) for ii, jj in zip(fem_list, ana_list)]
error_max = [np.max(ii) for ii in error_list]
print 'Error max for rad, tan, mix: ', error_max

# create color lists for plotting of potentials and error
vmax = .02
vmin = -vmax

# plotting
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
ax2.set_title('Analytical', y = 1.1)
ax3.set_title('FEM', y = 1.1)
ax4.set_title('Error', y = 1.1)


for ax, P1 in zip([ax1, ax5, ax9], P1s):
    color_list = ['#ffcc99', '#e6e6e6', '#83caff', '#999999']
    for idx in range(4):  # Hack - circles are overlaid on each other
        c = plt.Circle((0, 0), radii[-(idx + 1)], color=color_list[idx])
        ax.add_artist(c)
        # c_outline = plt.Circle((0, 0), radii[-(idx + 1)], color='gray',
        #                        lw=.3, fill=False)
        # ax.add_patch(c_outline)


    ax.arrow(-10, -10, 0, 1.8, head_width=.4, head_length=.4, fc='k', ec='k')
    ax.arrow(-10, -10, 1.8, 0, head_width=.4, head_length=.4, fc='k', ec='k')
    ax.arrow(-10, -10, -0.78, -0.78, head_width=.4, head_length=.33, fc='k', ec='k')
    plt.axis('equal')

    ax.text(-11.5, -12, 'x', size = 6.)
    ax.text(-7, -10.38, 'y', size = 6.)
    ax.text(-10.28, -7, 'z', size = 6.)

    ax.set_xlim(-12., 12.)
    ax.set_ylim(-12., 12.)

    rz1 = np.array((0, 0, 7.8))
    arrow = np.sum(P1, axis=0)*10
    start_pos = rz1 - arrow

    ax.arrow(start_pos[1], start_pos[2],
             2 * arrow[1], 2 * arrow[2],
             fc='k',
             ec='k',
             width=0.10,
             head_width=0.25,
             length_includes_head=False)
    ax.plot(rz1[1], rz1[2], 'ro', ms=4)
    ax.axis('off')


def clr(phi, error=False):
    if error is False:
        return plt.cm.PRGn((phi - vmin) / (vmax - vmin))
    else:
        return plt.cm.Greys((phi - vmin_error) / (vmax_error - vmin_error))

# create color lists for plotting of potentials and error
vmax = .05
vmin = -vmax
vmax_error = .001
vmin_error = 0.

# clr = lambda phi: plt.cm.PRGn((phi - vmin) / (vmax - vmin))
# clr_error = lambda error: plt.cm.Greys((error - vmin_error) / (vmax_error - vmin_error))

# colors_4s_list = [clr(ana_list[idx]) for idx in range(3)]
# colors_fem_list = [clr(fem_list[idx]) for idx in range(3)]
# colors_error_list = [clr_error(error_list[idx]) for idx in range(3)]


# plot analytically calculated potentials
X = params.x_points.reshape(180, 180)
Y = params.y_points.reshape(180, 180)
Z = params.z_points.reshape(180, 180)

for ax, ii in zip([ax2, ax6, ax10], [0, 1, 2]):
    surf1 = ax.plot_surface(X, Y, Z,
                            rstride=10, cstride=10, linewidth=0,
                            facecolors=clr(ana_list[ii]),
                            antialiased=False)

for ax, ii in zip([ax3, ax7, ax11], [0, 1, 2]):
    surf1 = ax.plot_surface(X, Y, Z,
                            rstride=10, cstride=10, linewidth=0,
                            facecolors=clr(fem_list[ii]),
                            antialiased=False)

for ax, ii in zip([ax4, ax8, ax12], [0, 1, 2]):
    surf1 = ax.plot_surface(X, Y, Z,
                            rstride=10, cstride=10, linewidth=0,
                            facecolors=clr(error_list[ii], error=True),
                            antialiased=False)


# # plot FEM-computed potentials
# for ax, clrs in zip([ax3, ax7, ax11], colors_fem_list):
#     surf2 = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=clrs,
#                              linewidth=0, antialiased=False)
# # plot errors
# for ax, clrs in zip([ax4, ax8, ax12], colors_error_list):
#     surf3 = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors = clrs,#cmap=plt.cm.Blues,
#                    linewidth=0, antialiased=False)  #, norm=colors.LogNorm(vmin=1e-7, vmax=1e-2))

for ax in [ax2, ax3, ax6, ax7, ax10, ax11, ax4, ax8, ax12]:  # [ax2, ax3, ax4, ax6, ax7, ax8, ax10, ax11, ax12]:
    ax.set_aspect('equal')
    ax.axis('off')
    # # ax.auto_scale_xyz([-60000., 60000.], [-60000., 60000.], [-60000., 60000.])
    ax.set_xlim3d(-6.5, 6.5)
    ax.set_ylim3d(-6.5, 6.5)
    ax.set_zlim3d(0.3-6.5, 0.3+6.5)
    ax.view_init(10, 0)

# # colorbars
# cax1 = fig.add_axes([0.31, 0.07, 0.38, 0.01])
# m = plt.cm.ScalarMappable(cmap=plt.cm.PRGn)
# m.set_array(phi_4s_list[1])
# cbar1 = fig.colorbar(m, cax=cax1, format='%3.3f', extend = 'both', orientation='horizontal')
# cbar1.outline.set_visible(False)
# ticks = np.linspace(-0.8, 0.8, 9)
# plt.xticks(ticks, [str(tick) for tick in ticks], rotation = 40)
# # multiply both P and phi's with 1000 and we get the same result, but in mV
# cbar1.set_label('Potential (mV)', labelpad=5.2)

# cax2 = fig.add_axes([0.79, 0.07, 0.15, 0.01])
# m = plt.cm.ScalarMappable(cmap=plt.cm.Greys)
# ticks2 = np.linspace(0.001, 0.005, 9)
# m.set_array(error_list[0])
# cbar2 = fig.colorbar(m, cax=cax2, format='%3.6f', extend='max', orientation='horizontal')
# cbar2.outline.set_visible(False)
# cax2.set_xticks(ticks2)
# cax2.set_xticklabels(['0.001', '', '0.002', '', '0.003', '', '0.004', '', '0.005'], rotation = 40)

# # multiply both P and phi's with 1000 and we get the same result, but in mV
# cbar2.set_label('Potential (mV)', labelpad=.1)

# fig.subplots_adjust(wspace = 0., hspace=-0.01, right = 0.98, left = 0.02)

# # label figures
# fig.text(0.03, .87, 'A',
#         horizontalalignment='center',
#         verticalalignment='center',
#         fontweight='demibold',
#         fontsize=12)
# fig.text(0.27, .87, 'B',
#         horizontalalignment='center',
#         verticalalignment='center',
#         fontweight='demibold',
#         fontsize=12)
# fig.text(0.51, .87, 'C',
#         horizontalalignment='center',
#         verticalalignment='center',
#         fontweight='demibold',
#         fontsize=12)
# fig.text(0.75, .87, 'D',
#         horizontalalignment='center',
#         verticalalignment='center',
#         fontweight='demibold',
#         fontsize=12)
# fig.text(0.03, .61, 'E',
#         horizontalalignment='center',
#         verticalalignment='center',
#         fontweight='demibold',
#         fontsize=12)
# fig.text(0.27, .61, 'F',
#         horizontalalignment='center',
#         verticalalignment='center',
#         fontweight='demibold',
#         fontsize=12)
# fig.text(0.51, .61, 'G',
#         horizontalalignment='center',
#         verticalalignment='center',
#         fontweight='demibold',
#         fontsize=12)
# fig.text(0.75, .61, 'H',
#         horizontalalignment='center',
#         verticalalignment='center',
#         fontweight='demibold',
#         fontsize=12)
# fig.text(0.03, .34, 'I',
#         horizontalalignment='center',
#         verticalalignment='center',
#         fontweight='demibold',
#         fontsize=12)
# fig.text(0.27, .34, 'J',
#         horizontalalignment='center',
#         verticalalignment='center',
#         fontweight='demibold',
#         fontsize=12)
# fig.text(0.51, .34, 'K',
#         horizontalalignment='center',
#         verticalalignment='center',
#         fontweight='demibold',
#         fontsize=12)
# fig.text(0.75, .34, 'L',
#         horizontalalignment='center',
#         verticalalignment='center',
#         fontweight='demibold',
#         fontsize=12)


# fig.set_size_inches(9., 6.)

# plt.savefig('./results/eeg_fig200000_new.pdf', dpi=600., bbox_inches='tight')
plt.show()
