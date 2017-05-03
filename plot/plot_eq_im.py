import sys
import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

class rayData:

    X = None
    xa = None
    ua = None
    xb = None
    ub = None

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

        self.N = N1*N2

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

        Xc = np.reshape(self.X[:,0], (self.N1, self.N2))
        Yc = np.reshape(self.X[:,0], (self.N1, self.N2))
        return Xc, Yc
    
    def _imageCartCoords_elps(self):

        R = np.reshape(self.X[:,0], (self.N1, self.N2))
        PHI = np.reshape(self.X[:,1], (self.N1, self.N2))

        Xc = R * np.cos(PHI)
        Yc = R * np.sin(PHI) * math.cos(self.i)
        return Xc, Yc





def plot_eq_i(filename, id):

    rays = rayData(filename)

    Inu = np.empty(rays.N)

    for i, 
        Inu[i] = 

    X, Y = rays.imageCartCoords()



if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("usage: $ python plot_eq_im.py <ray_files...>")

    for i, fname in enumerate(sys.argv[1:]):
        plot_eq_im(fname, i)
