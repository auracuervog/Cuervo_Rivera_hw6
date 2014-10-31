import matplotlib.pyplot as plt
import numpy as np
import sys
from mpl_toolkits.mplot3d import Axes3D

filename = sys.argv[1]
datos = np.loadtxt(filename)
x = datos[:,0]
y = datos[:,1]
z = datos[:,2]

fig = plt.figure()
ax = plt.axes()
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("Trayectoria xy")
plt.plot(x,y,'b')

plt.savefig(filename + '_2d'+ '.pdf',format = 'pdf', transparent=True)
plt.close()
