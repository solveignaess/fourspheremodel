from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import gridspec
from mpl_toolkits.axes_grid.inset_locator import inset_axes
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

abs_max_val = [np.max(np.abs(ii)) for ii in ana_list]
abs_max_range = [0.3 * ii for ii in abs_max_val]
print 'Error max for rad, tan, mix: ', error_max
print 'Abs max for rad, tan, min', [np.max(np.abs(ii)) for ii in ana_list]


def set_axis(ax, letter):
    ax.text(0.05,
            0.9,
            letter,
            fontsize=12,
            weight='bold',
            transform=ax.transAxes)
    return ax


def set_axis_3d(ax, letter):
    ax.text(0.05,
            1.025,
            100.025,
            letter,
            fontsize=12,
            weight='bold',
            transform=ax.transAxes)
    return ax


def draw_diagram(ax):
    color_list = ['#ffcc99', '#e6e6e6', '#83caff', '#999999']
    for idx in range(4):  # Hack - circles are overlaid on each other
        c = plt.Circle((0, 0), radii[-(idx + 1)], color=color_list[idx])
        ax.add_artist(c)
    ax.arrow(-10, -10, 0, 1.8, head_width=.4, head_length=.4, fc='k', ec='k')
    ax.arrow(-10, -10, 1.8, 0, head_width=.4, head_length=.4, fc='k', ec='k')
    ax.arrow(-10, -10, -0.78, -0.78, head_width=.4, head_length=.33, fc='k', ec='k')
    ax.text(-11.5, -12, 'x', size=6.)
    ax.text(-7, -10.38, 'y', size=6.)
    ax.text(-10.28, -7, 'z', size=6.)
    return ax


def adjust_3d_axis(ax):
    ax.set_xlim3d(-6.5, 6.5)
    ax.set_ylim3d(-6.5, 6.5)
    ax.set_zlim3d(0.3 - 6.5, 0.3 + 6.5)
    ax.view_init(10, 0)
    ax.set_aspect('equal')
    ax.axis('off')
    return ax


z_steps = 4
height_ratios = [1 for i in range(z_steps - 1)]
height_ratios.append(0.07)  # height of colorbar
fig = plt.figure()
gs = gridspec.GridSpec(z_steps, 4, height_ratios=height_ratios)

for ii, P1, letter in zip([0, 1, 2], P1s, ['A', 'E', 'I']):
    ax = plt.subplot(gs[ii, 0])
    draw_diagram(ax)
    rz1 = np.array((0, 0, 7.8))
    arrow = np.sum(P1, axis=0) * 10
    start_pos = rz1 - arrow
    ax.arrow(start_pos[1], start_pos[2],
             2 * arrow[1], 2 * arrow[2],
             fc='k',
             ec='k',
             width=0.10,
             head_width=0.25,
             length_includes_head=False)
    ax.plot(rz1[1], rz1[2], 'ro', ms=4)
    ax.set_xlim(-15., 15.)
    ax.set_ylim(-15., 15.)
    ax.axis('off')
    set_axis(ax, letter)

# plot analytically calculated potentials
X = params.x_points.reshape(180, 180)
Y = params.y_points.reshape(180, 180)
Z = params.z_points.reshape(180, 180)

rstride = 10
cstride = 10

title_texts = ['Analytical', 'FEM', 'Error']


def plot_phi(idx_val, letters, phi, error=None):
    for ii, letter in zip([0, 1, 2], letters):
        ax = plt.subplot(gs[ii, idx_val], projection='3d')
        if error is None:
            vmax = 0.5  # abs_max_range[ii]
            vmin = -1 * vmax
            clrs = plt.cm.PRGn((phi[ii] - vmin) / (vmax - vmin))
        else:
            vmax = 0.002
            vmin = 0.0
            clrs = plt.cm.Greys((phi[ii] - vmin) / (vmax - vmin))
        ax.plot_surface(X, Y, Z,
                        rstride=rstride, cstride=cstride, linewidth=0,
                        facecolors=clrs,
                        antialiased=False)
        if ii == 0:
            ax.set_title(title_texts[idx_val - 1], y=1.1)

        set_axis_3d(ax, letter)
        adjust_3d_axis(ax)


plot_phi(1, ['B', 'F', 'J'], ana_list)
plot_phi(2, ['C', 'G', 'K'], fem_list)
plot_phi(3, ['D', 'H', 'L'], error_list, error=True)

plt.show()
