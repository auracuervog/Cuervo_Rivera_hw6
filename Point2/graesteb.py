import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import sys
import matplotlib
filename = sys.argv[1]
datos = np.loadtxt(filename)
x=datos[:,0]
y=datos[:,1]
z=datos[:,2]


fig=plt.figure()
ax = fig.add_subplot(111, projection='3d')
o=fig.gca(projection='3d')
o.scatter(x,y,z)
o.set_xlabel('x')
o.set_ylabel('y')
o.set_zlabel('z')
plt.savefig(filename + '_3d'+ '.pdf',format = 'pdf', transparent=True)
plt.show() 
