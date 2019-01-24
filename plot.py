import numpy as np
import matplotlib.pyplot as plt
from sys import exit

rho = np.loadtxt("rho.dat")
bins = np.loadtxt("rho_bins.dat")

fx = np.loadtxt("fx.dat")
fy = np.loadtxt("fy.dat")
fz = np.loadtxt("fz.dat")

px = np.loadtxt("px.dat")
py = np.loadtxt("py.dat")
pz = np.loadtxt("pz.dat")

plt.subplot(2,3,1)
plt.imshow(rho,interpolation='none')

plt.subplot(2,3,2)
plt.imshow(fx,interpolation='none')

plt.subplot(2,3,3)
plt.imshow(fy,interpolation='none')

plt.subplot(2,3,5)
plt.imshow(px,interpolation='none')

plt.subplot(2,3,6)
plt.imshow(py,interpolation='none')

plt.show()


