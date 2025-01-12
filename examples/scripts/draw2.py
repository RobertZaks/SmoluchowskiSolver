#!/bin/python3
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

res = open('../graph.csv', "r")
tmp = res.readline().strip().split(',')
N = int(tmp[0])
M = int(tmp[1])
H, T, t1, t2 = [float(x) for x in tmp[2:]]
res.close()
res = np.loadtxt('graph.csv', delimiter=',', skiprows=1)

tau = T / M
m1 = math.floor(t1 / tau)
m2 = math.ceil(t2 / tau)
#print(N, M, H, T, t1, t2, m1, m2)
n = 8
whf = range(0, (m2 - m1 + 1) * (N + 1))
wh = np.linspace(0, H, N + 1, endpoint=True)
colors = cm.rainbow(np.linspace(0, 1, n + 1))

fig, ax = plt.subplots()
k = 2

m10 = math.floor(0.12 / tau)
m20 = math.floor(0.2 / tau)
m30 = math.floor(0.28 / tau)

ax.plot(res[(m10 - m1) * (N + 1):(m10 - m1 + 1) * (N + 1), 0], color=colors[0], label='v_0')
ax.plot(res[(m10 - m1) * (N + 1):(m10 - m1 + 1) * (N + 1), 1], color=colors[3], label='v_res(t=0.12)')
ax.plot(res[(m20 - m1) * (N + 1):(m20 - m1 + 1) * (N + 1), 1], color=colors[1], label='v_res(t=0.2)')
ax.plot(res[(m30 - m1) * (N + 1):(m30 - m1 + 1) * (N + 1), 1], color=colors[2], label='v_res(t=0.28)')
#ax.plot(res[(m10 - m1) * (N + 1):(m10 - m1 + 1) * (N + 1), 2], color=colors[8], linestyle='--', label='v_sol')
ax.plot(res[(m10 - m1) * (N + 1):(m10 - m1 + 1) * (N + 1), 2], color=colors[8], linestyle='--', label='v_sol(t=0.12)')
ax.plot(res[(m20 - m1) * (N + 1):(m20 - m1 + 1) * (N + 1), 2], color=colors[8], linestyle='--', label='v_sol(t=0.2)')
ax.plot(res[(m30 - m1) * (N + 1):(m30 - m1 + 1) * (N + 1), 2], color=colors[8], linestyle='--', label='v_sol(t=0.28)')

ax.legend(fontsize=10)
plt.xlabel('x', fontsize=10)
fig.savefig('resgraph.pdf', format='pdf', transparent=True)
