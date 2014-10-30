import matplotlib.pyplot as plt
import numpy as np
import sys

x0 = sys.argv[1]
filename = "presa%s.dat" %x0
datos = np.loadtxt(filename)
t=datos[:,0]
x=datos[:,1]
y=datos[:,2]

fig = plt.figure()
ax = plt.axes()
ax.set_xlabel("t")
ax.set_ylabel("x")
ax.set_title("Presas")
plt.plot(t,x,'r',label="Presas")
ax.legend()
plt.plot(t,y,'b',label="Cazadores")
ax.legend()

filename = 'presas_x'
plt.savefig(filename + x0 +'.pdf',format = 'pdf', transparent=True)
plt.close()


