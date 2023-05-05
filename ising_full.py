import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("Ising_full.dat")

time = data[:, 0]
energy = data[:, 1]
magnetization = data[:, 2]

plt.plot(time, magnetization, linewidth=0.2, color='black', label='<M>/N')
plt.plot(time, energy, linewidth=0.2, color='lime', label='<E>/N')
plt.legend()
plt.ylabel('var_E/N and var_M/N')
plt.xlabel('T')
plt.show()
