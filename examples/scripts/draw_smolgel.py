#!/bin/python3
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

res = np.loadtxt('../smolgel_conc.csv', delimiter='\t')
fig, ax = plt.subplots()
ax.plot(res[:, 0], res[:, 1], color='red', label='c')
ax.legend(fontsize=10)
plt.grid()
plt.xlabel('x', fontsize=10)
plt.ylabel('x * c', fontsize=10)
fig.savefig('smolgel_conc.pdf', format='pdf')

res = np.loadtxt('../smolgel_mass.csv', delimiter='\t')
fig, ax = plt.subplots()
ax.plot(res[:, 0], res[:, 1], color='green', label='mass')
N = len(res[:, 0])
T = 2
k = int(N / T)
wt = np.linspace(0, T, N, endpoint=True)
#ax.plot(wt[k:], list(map(lambda z: 1+0, wt[k:])), color='red', label='1')
ax.legend(fontsize=10)
plt.grid()
plt.xlabel('t', fontsize=10)
fig.savefig('smolgel_mass.pdf', format='pdf')
