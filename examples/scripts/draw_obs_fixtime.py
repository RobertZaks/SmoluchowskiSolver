#!/bin/python3
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

res = np.loadtxt('../graph_obs_fixtime.csv', delimiter=',', skiprows=1)
fig, ax = plt.subplots()
ax.plot(res[:, 0], color='red', label='c_obs_s')
ax.plot(res[:, 1], color='green', linestyle='--', label='c_obs_r')
ax.legend(fontsize=10)
plt.grid(axis='x')
plt.xlabel('x', fontsize=10)
fig.savefig('resgraph_obs_fixtime.pdf', format='pdf', transparent=True)
