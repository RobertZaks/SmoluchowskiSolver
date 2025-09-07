#!/bin/python3
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

res = np.loadtxt('../graph_obs.csv', delimiter=',')
res2 = np.loadtxt('../graph_obs_fixtime.csv', delimiter=',')
fig, ax = plt.subplots()
N = len(res[:,0]) - 1
H = 6.0
whf = H / N  * np.arange(N + 1)
#ax.plot(whf, res[:, 0], color='0.25', label='c=c_obs_s')
#ax.plot(whf, res[:, 1], color='0', linestyle='--', label='c=c_obs_r')
ax.plot(whf, res[:, 0], color='red', label='c=c_obs_s')
ax.plot(whf, res[:, 1], color='green', linestyle='--', label='c=c_obs_r')
ax.plot(whf, res2[:, 1], color='blue', linestyle='--', linewidth=2, label='c=c_obs_r_fixtime')
ax.legend(fontsize=10)
plt.grid()
#plt.yscale('log')
plt.xlabel('x', fontsize=10)
plt.ylabel('x * c', fontsize=10)
#fig.savefig('resgraph_obs.eps', format='eps')
fig.savefig('resgraph_obs.pdf', format='pdf', transparent=True)
