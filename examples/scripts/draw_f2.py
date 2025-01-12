#!/bin/python3
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

fig, ax = plt.subplots()
colors = ['blue', 'red', 'green', 'yellow']
for i, alpha in enumerate(['05', '01', '005', '001']):
    res = np.loadtxt('../functional' + alpha + '.csv', delimiter=',')
    ax.plot(res, color=colors[i], label='alpha=0.' + alpha)
    ymin = res.min()
    xmin = np.where(res == ymin)[0].min()
    ax.plot(xmin, ymin, color=colors[i], marker='o')
    ax.axhline(ymin, linestyle='--', color=colors[i])
ax.legend(fontsize=10)
plt.grid()
#plt.yscale('log')
plt.xlabel('num_iter', fontsize=10)
plt.ylabel('J(v)', fontsize=10)
fig.savefig('resgraph_functional.pdf', format='pdf', transparent=True)
