import math
import sys
import matplotlib.pyplot as plt
import numpy as np

def plot_t(file, ax, *args, **kwargs):

    t, x, v = np.loadtxt(file, unpack=True, usecols=[0,1,2])

    ax.plot(t, x, *args, **kwargs)


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print "\nPlease provide a data file\n"
        sys.exit()

    color = ['k', 'b', 'g', 'r']

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for i,f in enumerate(sys.argv[1:]):
        plot_t(f, ax, ls='-', marker='.', color=color[i], 
                label="{0}".format(f))
    xlim = ax.get_xlim()
    T = np.linspace(xlim[0], xlim[1], num=10000)
    ax.plot(T, np.cos(T), color='grey')
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$x$")
    #ax.set_aspect('equal')
    ax.legend()

    plt.show()
