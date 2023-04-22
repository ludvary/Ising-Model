import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("Ising_full.dat")

time = data[:, 0]
magnetization = data[:, 1]
energy = data[:, 2]

plt.plot(time, energy, linewidth=0.2, color='black', label='<M>/N')
plt.plot(time, magnetization, linewidth=0.2, color='lime', label='<E>/N')
plt.legend()
plt.ylabel('avg_E/N and avg_M/N')
plt.xlabel('T')
plt.show()