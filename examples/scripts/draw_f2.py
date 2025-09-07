#!/bin/python3
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

plt.rcParams['text.usetex'] = True

fig, ax = plt.subplots()
colors = cm.rainbow(np.linspace(0, 1, 6))
#colors = ['blue', 'red', 'green', 'yellow']
base = '0.0'
num = 50
#for i in range(1, 10):
    #alpha = base + str(i)
for i, alpha in enumerate(['0.001', '0.005', '0.01', '0.02', '0.03', '0.05']):
    res = np.loadtxt('../res/functional' + alpha + '.csv', delimiter=',')
    res = res[:num]
    ax.plot(range(num), res, color=colors[i], label=r'$\alpha$' + '=' + alpha)
    ymin = res.min()
    xmin = np.where(res == ymin)[0].min()
    #ax.plot(xmin, ymin, color=colors[i], marker='o')
    ax.axhline(ymin, linestyle='--', color=colors[i])
ax.legend(fontsize=10)
plt.grid()
#plt.yscale('log')
plt.xlabel('num_iter', fontsize=10)
plt.ylabel(r'$J_{\alpha}(v)$', fontsize=10)
fig.savefig('resgraph_functional.pdf', format='pdf', transparent=True)
