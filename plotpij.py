import numpy as np
import matplotlib.pyplot as plt

fil = open('pij.dat')
vals = fil.readlines()
fil.close()
y = np.zeros(len(vals))
x = np.zeros(len(vals))
for i in range(len(vals)):
    value = vals[i].split()
    x[i] = float(value[0])
    y[i] = float(value[1])

print((y[1]-y[0])/(x[1]-x[0]))
plt.plot(x,y)
plt.xlabel('$\\theta$ [deg]')
plt.ylabel('$p_E(\\theta)$')
plt.savefig('pij.pdf')
