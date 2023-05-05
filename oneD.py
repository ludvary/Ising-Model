from matplotlib import pyplot as plt
import numpy as np

data = np.loadtxt("oneD.dat")

temperature = data[:, 0]
energy = data[:, 1]
magnetization = data[:, 2]

plt.plot(temperature, energy, linewidth='0.1', label='var_E', color='lime')
plt.plot(temperature, magnetization, linewidth='0.1', label='var_M', color='black')
plt.legend()
plt.ylabel('var_E/ var_M')
plt.xlabel('T')
plt.show()

