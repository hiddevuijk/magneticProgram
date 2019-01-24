from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from sys import exit

r = np.loadtxt("final_config.dat")

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')

for i in range(r.shape[0]):
	ax.scatter(r[i,0],r[i,1],r[i,2])

ax.set_xlim([0,10])
ax.set_ylim([0,10])
ax.set_zlim([0,10])

plt.show()



