#!/bin/python3
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

fig, ax = plt.subplots()
colors = cm.rainbow(np.linspace(0, 1, 6 + 1))
res = np.loadtxt('../graph_obs300.csv', delimiter=',')
ax.plot(res[:, 0], color=colors[0], linewidth=5, label=f'c_obs_s')
for i, num_iter in enumerate(['300', '600', '900']):
    res = np.loadtxt('graph_obs' + num_iter + '.csv', delimiter=',')
    res2 = np.loadtxt('graph_obs_fixtime' + num_iter + '.csv', delimiter=',')
    ax.plot(res[:, 1], color=colors[2 * i + 1], linestyle='--',
            label=f'c_obs_r ({num_iter =})')
    ax.plot(res2[:, 1], color=colors[2 * i + 2], linestyle='--', linewidth=2,
            label=f'c_obs_r_fixtime ({num_iter =})')
ax.legend(fontsize=10)
plt.grid()
plt.title('a)')
plt.yscale('log')
plt.xlabel('x', fontsize=10)
fig.savefig('resgraph_obs.pdf', format='pdf', transparent=True)
