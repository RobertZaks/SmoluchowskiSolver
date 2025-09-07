#!/bin/python3
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

res = np.loadtxt('../graph_smoltest.csv', delimiter='\t')
fig, ax = plt.subplots()
ax.plot(res[:, 0], res[:, 1], color='red', label='res')
ax.plot(res[:, 0], res[:, 2], color='blue', linestyle='--', label='sol')
ax.legend(fontsize=10)
plt.yscale('log')
plt.grid()
plt.xlabel('x', fontsize=10)
fig.savefig('resgraph_smoltest.pdf', format='pdf')
