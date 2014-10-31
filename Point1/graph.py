import matplotlib.pyplot as plt
import numpy as np
import sys

x0 = sys.argv[1]
filename = "presa%s.txt" %x0
datos = np.loadtxt(filename)
t=datos[:,0]
x=datos[:,1]
y=datos[:,2]

fig = plt.figure()
ax = plt.axes()
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("Presas")
plt.plot(x,y,'b',label="Presas")
ax.legend()

filename = 'presas_x'
plt.savefig(filename + x0 +'.pdf',format = 'pdf', transparent=True)
plt.close()


