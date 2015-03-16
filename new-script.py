from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("arrow-file.txt")
print data[:,1]
sorted_data = np.sort(data[:,1])
plt.plot(sorted_data[::-1])
plt.show()
