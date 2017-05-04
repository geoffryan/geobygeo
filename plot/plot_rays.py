import math
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def plot_ray(filename, id):

    s = np.loadtxt(filename, skiprows=1, usecols=[0])
    x = np.loadtxt(filename, skiprows=1, usecols=[1,2,3,4])
    u = np.loadtxt(filename, skiprows=1, usecols=[5,6,7,8])

    sa = s.min()
    sb = s.min() + 10
    ind = (s>=sa)*(s<=sb)

    s = s[ind]
    x = x[ind,:]
    u = u[ind,:]

    fig, ax = plt.subplots(2,2, figsize=(12,9))
    ax[0,0].plot(s, x[:,0])
    ax[0,0].set_xlabel(r'$\lambda$')
    ax[0,0].set_ylabel(r'$x^0$')
    ax[0,1].plot(s, x[:,1])
    ax[0,1].set_xlabel(r'$\lambda$')
    ax[0,1].set_ylabel(r'$x^1$')
    ax[1,0].plot(s, x[:,2])
    ax[1,0].set_xlabel(r'$\lambda$')
    ax[1,0].set_ylabel(r'$x^2$')
    ax[1,1].plot(s, x[:,3])
    ax[1,1].set_xlabel(r'$\lambda$')
    ax[1,1].set_ylabel(r'$x^3$')

    fig.tight_layout()

    outname = ".".join(filename.split(".")[:-1]) + ".png"
    fig.savefig(outname)
    plt.close(fig)

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("usage: $ python plot_rays.txt <rayfiles...>")
        sys.exit()

    for id, filename in enumerate(sys.argv[1:]):
        plot_ray(filename, id)

