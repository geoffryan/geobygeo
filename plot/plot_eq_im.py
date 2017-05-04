import sys
import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

class rayData:

    xa = None
    ua = None
    xb = None
    ub = None
    xc = None

    def __init__(self, filename=None):
        self.metric = -1
        self.grid   = -1
        self.N1     = 0
        self.N2     = 0
        self.X1lim  = (0,0)
        self.X2lim  = (0,0)
        self.D      = -1
        self.i      = -1
        self.az     = -1
        self.N = 0

        if filename is not None:
            self.loadRayFile(filename)

    def loadRayFile(self, filename):
        f = open(filename, "r")
        f.readline()
        self.metric = int(f.readline().split()[-1])
        self.grid = int(f.readline().split()[-1])
        self.N1 = int(f.readline().split()[-1])
        self.N2 = int(f.readline().split()[-1])
        f.readline()
        f.readline()
        f.readline()
        f.readline()
        self.D = float(f.readline().split()[-1])
        self.i = float(f.readline().split()[-1])
        self.az = float(f.readline().split()[-1])
        f.close()

        self.N = self.N1*self.N2

        self.xc = np.loadtxt(filename, skiprows=12, usecols=[2,3])
        self.xb = np.loadtxt(filename, skiprows=12, usecols=[4,5,6,7])
        self.ub = np.loadtxt(filename, skiprows=12, usecols=[8,9,10,11])
        self.xa = np.loadtxt(filename, skiprows=12, usecols=[12,13,14,15])
        self.ua = np.loadtxt(filename, skiprows=12, usecols=[16,17,18,19])

    def imageCartCoords(self):

        if self.grid == 0:
            return self._imageCartCoords_cart()
        elif self.grid == 1:
            return self._imageCartCoords_elps()

        return None, None

    def _imageCartCoords_cart(self):

        Xc = np.reshape(self.xc[:,0], (self.N1, self.N2))
        Yc = np.reshape(self.xc[:,1], (self.N1, self.N2))
        return Xc, Yc
    
    def _imageCartCoords_elps(self):

        R = np.reshape(self.xc[:,0], (self.N1, self.N2))
        PHI = np.reshape(self.xc[:,1], (self.N1, self.N2))

        Xc = R * np.cos(PHI)
        Yc = R * np.sin(PHI) * math.cos(self.i)
        return Xc, Yc


def spiral(nu, x, u, metric):

    m = 2

    if metric == 0:
        r = math.sqrt(x[1]*x[1] + x[2]*x[2])
        phi = math.atan2(x[2], x[1])
        us = np.array([1.0, 0.0, 0.0, 0.0])
    elif metric == 1:
        r = x[1] * math.sin(x[2])
        phi = x[3]
        us = np.array([1.0, 0.0, 0.0, 0.0])
    elif metric == 2:
        r = x[1] * math.sin(x[2])
        phi = x[3]
        if x[1] > 0.0:
            us = np.array([1.0/math.sqrt(1-2/x[1]), 0.0, 0.0, 0.0])
        else:
            us = np.array([0.0, 0.0, 0.0, 0.0])
    elif metric == 3:
        r = x[1] * math.sin(x[2])
        phi = x[3]
        if x[1] > 0.0:
            us = np.array([1.0/math.sqrt(1-2/x[1]), 0.0, 0.0, 0.0])
        else:
            us = np.array([0.0, 0.0, 0.0, 0.0])
    
    PHI = m * (phi - 2*math.pi*r/(30.0))

    return math.cos(0.5*PHI)**4 * np.ones(nu.shape)

def plot_eq_im(filename, id, imfunc):

    rays = rayData(filename)

    Inu = []

    nu = np.array([1.0])

    for i in range(rays.N):
        Inu.append(imfunc(nu, rays.xa[i], rays.ua[i], rays.metric))

    Inu = np.array(Inu)
    #Inu.reshape((rays.N1, rays.N2, -1))

    X, Y = rays.imageCartCoords()

    tri = mpl.tri.Triangulation(X.flat, Y.flat)

    for j in range(len(nu)):
        fig, ax = plt.subplots(1,1, figsize=(8,6))
        ax.tricontourf(tri, Inu[:,j], 256, cmap=mpl.cm.inferno)
        #ax.plot(X[:,:-1].flat, Y[:,:-1].flat, ls='', marker='.', ms=10, color='green')
        ax.set_aspect('equal')
        fig.savefig("im_{0:03d}_{1:01d}.png".format(id, j), dpi=300)
        plt.close(fig)


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("usage: $ python plot_eq_im.py <ray_files...>")

    for i, fname in enumerate(sys.argv[1:]):
        plot_eq_im(fname, i, spiral)
