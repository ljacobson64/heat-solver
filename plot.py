#!/usr/bin/python

import numpy as np
import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib import rc
from mpl_toolkits.mplot3d import Axes3D

def main():
    to_file = True
    combined = False
    dpi = 300

    rc('font', **{'size': 10, 'family': 'serif', 'serif': ['Computer Modern Roman']})
    rc('text', usetex=True)

    # Load data from file
    data_fwd = np.loadtxt('T_fwd.txt')
    data_adj = np.loadtxt('T_adj.txt')

    # Setup grid
    Lx = 2
    Ly = 1
    nx = data_fwd.shape[1] - 1
    ny = data_fwd.shape[0] - 1
    x = np.linspace(0, Lx, nx + 1)
    y = np.linspace(0, Ly, ny + 1)
    xp, yp = np.meshgrid(x, y)

    Lmax = max(Lx, Ly)

    if combined: fig = plt.figure(figsize=(24, 6))

    for i in range(3):
        if combined: ax = fig.add_subplot(1, 3, i, projection='3d')
        else: fig = plt.figure()
        ax = fig.gca(projection='3d')
        if   i == 0: data = data_fwd
        elif i == 1: data = data_adj
        elif i == 2: data = data_fwd * data_adj
        ax.plot_surface(xp, yp, data, rstride = 1, cstride = 1,
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
            ax.set_zlabel(r'Temperature [K]')
            ax.set_title(r'Forward solution')
            if to_file: plt.savefig('forward.png', dpi=dpi); plt.close()
        elif i == 1:
            ax.set_zlabel(r'Adjoint temperature $\left[\frac{\textup{m}^3\textup{-K}}{\textup{W}}\right]$')
            ax.set_title(r'Adjoint solution')
            if to_file: plt.savefig('adjoint.png', dpi=dpi); plt.close()
        elif i == 2:
            ax.set_zlabel(r'Forward $\times$ adjoint')
            ax.set_title(r'Contributon')
            if to_file: plt.savefig('contributon.png', dpi=dpi); plt.close()

    if not to_file: plt.show()

if __name__ == '__main__':
    main()
