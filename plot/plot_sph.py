import math
import sys
import matplotlib.pyplot as plt
import numpy as np

def plot_xy(file, ax, *args, **kwargs):

    ta, t, r, th, phi = np.loadtxt(file, unpack=True, usecols=[0,1,2,3,4], 
                                    skiprows=1)

    x = r*np.sin(th)*np.cos(phi)
    y = r*np.sin(th)*np.sin(phi)

    ax.plot(x, y, *args, **kwargs)


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print "\nPlease provide a data file\n"
        sys.exit()

    color = ['k', 'b', 'g', 'r']

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for i,f in enumerate(sys.argv[1:]):
        plot_xy(f, ax, ls='-', marker='.', color=color[i], 
                label="{0}".format(f))
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")
    ax.set_xlim(0,1)
    ax.set_ylim(0,2)
    ax.set_aspect('equal')
    ax.legend()

    plt.show()
