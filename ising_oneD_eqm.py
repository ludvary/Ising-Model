import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt('oneD_eqm.dat')

t = data[:, 0]
E = data[:, 1]
M = data[:, 2]

plt.plot(t, E, linewidth='0.1', color='teal', label='E/N')
plt.plot(t, M, linewidth='0.1', color='purple', label='M/N')
plt.legend()
plt.ylabel('E/N and M/N')
plt.xlabel('T')
plt.show()
