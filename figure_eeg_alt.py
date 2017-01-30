from __future__ import division
import numpy as np
import math
from CalcPotential4Sphere import CalcPotential4Sphere
import matplotlib.pyplot as plt
from matplotlib import colors
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

plt.close('all')
file_eeg = np.load('./results/chaitanya_wm_sigma20_poly2_eeg.npz')

files = [file_eeg]  # [file_res10, file_res50, file_res100, file_res200] #  [file_res10, file_res30, file_res50, file_res70, file_res100]

for f in files:
    param_dict = f['params'].item()
    fem_phi = f['phi']
    # X = f['X']
    # Y = f['Y']
    # Z = f['Z']
    charge_pos = param_dict['charge_pos']
    charges = param_dict['charge']
    rz1 = (charge_pos[0] + charge_pos[1])/2
    # rz1 = np.array([0., 0., 10000.])
    P1 = np.array([charge_pos[0]*charges[0] + charge_pos[1]*charges[1]])

    r = param_dict['radii'][3]
    radii = param_dict['radii']
    sigmas = param_dict['sigmas']
    # sigmas[1] += 1e-14

    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    X = (r-1e2) * np.outer(np.cos(u), np.sin(v))
    Y = (r-1e2) * np.outer(np.sin(u), np.sin(v))
    Z = (r-1e2) * np.outer(np.ones(np.size(u)), np.cos(v))

    RE = np.zeros(X.shape)
    phi_4s = np.zeros(X.shape)
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            sphere_mod = CalcPotential4Sphere(radii, sigmas, np.array([X[i,j], Y[i,j], Z[i,j]]), rz1)
            pot = sphere_mod.calc_potential(P1)
            phi_4s[i, j] = np.nan_to_num(pot)

    RE = np.abs(fem_phi - phi_4s) / np.abs(phi_4s)
    xs = np.linspace(-r, r, 201)
    zs = np.linspace(-r, r, 201)
    elec_pos = []
    for zpos in zs:
        for xpos in xs:
            elec_pos.append([xpos, 0., zpos])
    elec_pos = np.array(elec_pos)
    Xm, Zm = np.meshgrid(xs, zs)
    Ym = np.zeros(X.shape)

    subdomain_markers = []
    for pos in elec_pos:
        r = np.sqrt(pos[0]**2 + pos[1]**2 + pos[2]**2)
        #print r
        if r < radii[0]:  #  or (r-79000.) < 1e-16:
            subdomain_markers.append(0.1)
        elif r < radii[1]: #  or (r-80000.) < 1e-16:
            subdomain_markers.append(0.2)
        elif r < radii[2]: #  or (r-90000.) < 1e-16:
            subdomain_markers.append(0.3)
        elif r < radii[3]:  # or (r-100000.) < 1e-16:
            subdomain_markers.append(0.4)
        else:
            subdomain_markers.append(np.nan)
    subd_markers = np.array(subdomain_markers).reshape(Xm.shape)

    #plt.close('all')
    fig = plt.figure(figsize=(12, 4))

    ax1 = plt.subplot(1, 4, 1)
    im = ax1.imshow(subd_markers.reshape(len(xs), len(zs)).T,
                    origin='lower', interpolation='nearest', extent = [-90000., 90000., -90000., 90000.])
    cmap = colors.ListedColormap([ '#FEC0C2', '#FFA93A', '#B3B0FF', '#FFC072'])#[colorbrewer['green'], colorbrewer['lightgreen']])
    im.set_cmap(cmap)
    ax1.add_patch(plt.Circle((0, 0), radius = radii[0], color = 'gray', fill=False, lw = .3))
    ax1.add_patch(plt.Circle((0, 0), radius = radii[1], color = 'gray', fill=False, lw = .3))
    ax1.add_patch(plt.Circle((0, 0), radius = radii[2], color = 'gray', fill=False, lw = .3))
    ax1.add_patch(plt.Circle((0,0), radius = radii[3], color = 'gray', fill=False, lw = .5))
    arrow = np.sum(P1, axis = 0)
    ax1.arrow(charge_pos[0][0], charge_pos[0][2],
              4*arrow[0], 4*arrow[2],
              fc='k',  #colorbrewer['green'],
              ec='k',  # colorbrewer['green'],
              width = 200,
              length_includes_head=False)
    ax1.axis('off')

    vmax = np.max(np.abs(phi_4s))
    vmin = -vmax
    xx, yy = np.meshgrid(u, v)
    levels = np.arange(vmin, vmax, (vmax - vmin) / 18.)
    plt.subplot(1, 4, 2)
    img = plt.contourf(xx, yy, fem_phi, levels=levels, cmap=plt.cm.PRGn)
    plt.xlabel('Azimulthal angle (degrees)')
    plt.ylabel('Polar angle (degrees)')
    plt.xticks(np.linspace(0, 2 * np.pi, 6), np.linspace(0, 360, 6))
    plt.yticks(np.linspace(0, np.pi, 6), np.linspace(0, 180, 6))
    plt.title('FEM')
    plt.colorbar(orientation='horizontal')

    plt.subplot(1, 4, 3)
    img = plt.contourf(xx, yy, phi_4s, levels=levels, cmap=plt.cm.PRGn)
    plt.xlabel('Azimulthal angle (degrees)')
    plt.ylabel('Polar angle (degrees)')
    plt.xticks(np.linspace(0, 2 * np.pi, 6), np.linspace(0, 360, 6))
    plt.yticks(np.linspace(0, np.pi, 6), np.linspace(0, 180, 6))
    plt.title('Analytical')
    plt.colorbar(orientation='horizontal')

    plt.subplot(1, 4, 4)
    img = plt.contourf(xx, yy, RE, 18, cmap=plt.cm.Blues)
    plt.xlabel('Azimulthal angle (degrees)')
    plt.ylabel('Polar angle (degrees)')
    plt.xticks(np.linspace(0, 2 * np.pi, 6), np.linspace(0, 360, 6))
    plt.yticks(np.linspace(0, np.pi, 6), np.linspace(0, 180, 6))
    plt.title('Relative Error')
    plt.colorbar(orientation='horizontal')

    #plt.tight_layout()

    #plt.colorbar()
    plt.show()
