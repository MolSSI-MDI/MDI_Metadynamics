import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('work/bias.dat', dtype=float, comments='#', usecols=(0,1))
print(len(data[:, 1]))

plt.plot(data[:,1])
plt.show()
