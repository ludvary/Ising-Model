import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("lattice_eqm.dat")

time = data[:, 0]
magnetization = data[:, 1]
energy = data[:, 2]

plt.plot(time, magnetization, linewidth=0.9, color='purple', label='M/N')
plt.plot(time, energy, linewidth=0.9, color='teal', label='E/N')
plt.legend()
plt.ylabel('E/N and M/N')
plt.xlabel('T')
plt.show()

