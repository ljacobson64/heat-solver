#!/usr/bin/python

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors, cm, rc
from mpl_toolkits.mplot3d import Axes3D

def main():
    to_file = True
    combined = False
    dpi = 300

    rc('font', **{'size': 10, 'family': 'serif', 'serif': ['Computer Modern Roman']})
    rc('text', usetex=True)

    # Load data from file
    data = []
    data.append(np.loadtxt('output/k.txt'))
    data.append(np.loadtxt('output/Q_fwd.txt'))
    data.append(np.loadtxt('output/Q_adj.txt'))
    data.append(np.loadtxt('output/T_fwd.txt'))
    data.append(np.loadtxt('output/T_adj.txt'))
    data.append(data[3] * data[4])
    num_plots = len(data)

    # Setup grid
    Lx = 2
    Ly = 1
    nx = data[3].shape[1] - 1
    ny = data[3].shape[0] - 1
    dx = Lx / nx
    dy = Ly / ny
    x_int = np.linspace(0.5 * dx, Lx - 0.5 * dx, nx)
    y_int = np.linspace(0.5 * dy, Ly - 0.5 * dy, ny)
    x_nod = np.linspace(0, Lx, nx + 1)
    y_nod = np.linspace(0, Ly, ny + 1)
    xp_int, yp_int = np.meshgrid(x_int, y_int)
    xp_nod, yp_nod = np.meshgrid(x_nod, y_nod)

    Lmax = max(Lx, Ly)

    if not os.path.isdir('img'): os.makedirs('img')

    if combined: fig = plt.figure(figsize=(8*num_plots, 6))

    for i in range(num_plots):
        if combined: ax = fig.add_subplot(1, num_plots, i, projection='3d')
        else: fig = plt.figure()
        ax = fig.gca(projection='3d')
        if   i in [0, 1, 2]: xp = xp_int; yp = yp_int
        elif i in [3, 4, 5]: xp = xp_nod; yp = yp_nod
        ax.plot_surface(xp, yp, data[i], rstride = 1, cstride = 1,
                        linewidth = 0, cmap = cm.get_cmap('Spectral_r'))
        ax.set_xlim([-0.5*(Lmax - Lx), Lx + 0.5*(Lmax - Lx)])
        ax.set_ylim([-0.5*(Lmax - Ly), Ly + 0.5*(Lmax - Ly)])
        xlims = ax.get_xlim()
        ylims = ax.get_ylim()
        zlims = ax.get_zlim()
        ax.set_zlim([0, zlims[1]])
        ax.ticklabel_format(axis='z', style='sci', scilimits=(0, 0))
        ax.zaxis.get_offset_text().set_visible(False)
        e_str = r'$10^{%u}$' % int(np.floor(np.log10(zlims[1])))
        ax.text(1.05 * xlims[1], 1.05 * ylims[1], 1.05 * zlims[1], e_str)
        ax.set_xlabel(r'x [m]')
        ax.set_ylabel(r'y [m]')
        if i == 0:
            ax.set_title(r'Thermal conductivity')
            ax.set_zlabel(r'Thermal conductivity [W/m-K]')
            if to_file: plt.savefig('img/k.png', dpi=dpi); plt.close()
        elif i == 1:
            ax.set_title(r'Heat source')
            ax.set_zlabel(r'Volumetric heat source [W/m^3]')
            if to_file: plt.savefig('img/forward_source.png', dpi=dpi); plt.close()
        elif i == 2:
            ax.set_title(r'Adjoint source')
            ax.set_zlabel(r'Volumetric adjoint source [-]')
            if to_file: plt.savefig('img/adjoint_source.png', dpi=dpi); plt.close()
        elif i == 3:
            ax.set_title(r'Forward solution')
            ax.set_zlabel(r'Temperature [K]')
            if to_file: plt.savefig('img/forward_soln.png', dpi=dpi); plt.close()
        elif i == 4:
            ax.set_title(r'Adjoint solution')
            ax.set_zlabel(r'Adjoint temperature $\left[\frac{\textup{m}^3\textup{-K}}{\textup{W}}\right]$')
            if to_file: plt.savefig('img/adjoint_soln.png', dpi=dpi); plt.close()
        elif i == 5:
            ax.set_title(r'Contributon')
            ax.set_zlabel(r'Forward $\times$ adjoint')
            if to_file: plt.savefig('img/contributon.png', dpi=dpi); plt.close()

    if not to_file: plt.show()

if __name__ == '__main__':
    main()
