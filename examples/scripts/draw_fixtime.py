#!/bin/python3
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

colors = cm.rainbow(np.linspace(0, 1, 9))
res = np.loadtxt('../graph_fixtime.csv', delimiter=',', skiprows=1)
fig, ax = plt.subplots()
ax.plot(res[:, 0], color=colors[0], label='v_0')
ax.plot(res[:, 1], color=colors[1], label='v_res_fixtime')
ax.plot(res[:, 2], color=colors[8], linestyle='--', label='v_sol')
ax.legend(fontsize=10)
plt.xlabel('x', fontsize=10)
fig.savefig('resgraph_fixtime.pdf', format='pdf', transparent=True)
