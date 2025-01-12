#!/bin/python3
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

res = np.loadtxt('../functional.csv', delimiter=',')
fig, ax = plt.subplots()

ymin = res.min()
xmin = np.where(res == ymin)
ax.plot(res[:], color='blue', label='J(v)')
ax.plot(xmin, ymin, color='yellow', marker='o')
ax.axhline(ymin, color='red')
ax.legend(fontsize=10)
plt.grid()
plt.yscale('log')
plt.xlabel('num_iter', fontsize=10)
fig.savefig('resgraph_functional.pdf', format='pdf', transparent=True)
