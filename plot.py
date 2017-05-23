#!/usr/bin/python

import numpy as np
import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def main():
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

    to_file = True
    combined = False

    if combined: fig = plt.figure(figsize=(24, 6))

    if combined: ax1 = fig.add_subplot(1, 3, 1, projection='3d')
    else: fig = plt.figure()
    ax1 = fig.gca(projection='3d')
    ax1.plot_surface(xp, yp, data_fwd, rstride = 1, cstride = 1,
                     linewidth = 0, cmap = cm.get_cmap('Spectral_r'))
    ax1.set_xlabel('x [m]')
    ax1.set_ylabel('y [m]')
    ax1.set_zlabel('Temperature [K]')
    ax1.set_title('Forward solution')
    if to_file: plt.savefig('forward.png', dpi=300); plt.close()

    if combined: ax1 = fig.add_subplot(1, 3, 2, projection='3d')
    else: fig = plt.figure()
    ax2 = fig.gca(projection='3d')
    ax2.plot_surface(xp, yp, data_adj, rstride = 1, cstride = 1,
                     linewidth = 0, cmap = cm.get_cmap('Spectral_r'))
    ax2.set_xlabel('x [m]')
    ax2.set_ylabel('y [m]')
    ax2.set_zlabel('Adjoint temperature [m^3-K/W]')
    ax2.set_title('Adjoint solution')
    if to_file: plt.savefig('adjoint.png', dpi=300); plt.close()

    if combined: ax1 = fig.add_subplot(1, 3, 3, projection='3d')
    else: fig = plt.figure()
    ax3 = fig.gca(projection='3d')
    ax3.plot_surface(xp, yp, data_fwd*data_adj, rstride = 1, cstride = 1,
                     linewidth = 0, cmap = cm.get_cmap('Spectral_r'))
    ax3.set_xlabel('x [m]')
    ax3.set_ylabel('y [m]')
    ax3.set_zlabel('Forward * adjoint')
    ax3.set_title('Contributon')
    if to_file: plt.savefig('contributon.png', dpi=300); plt.close()

    if not to_file: plt.show()

if __name__ == '__main__':
    main()
