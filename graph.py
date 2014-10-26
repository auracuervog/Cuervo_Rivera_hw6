import matplotlib.pyplot as plt
import numpy as np
datos = np.loadtxt('presa.dat')
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

filename = 'Presas'
plt.savefig(filename + '.pdf',format = 'pdf', transparent=True)
plt.close()


