import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("spin_flips_in_mass.dat")

time = data[:, 0]
magnetization = data[:, 1]

plt.plot(time, magnetization, linewidth=0.1, color='black' , label='<M/N')
plt.legend()
plt.ylabel('M/N')
plt.xlabel('MC steps')
plt.show()