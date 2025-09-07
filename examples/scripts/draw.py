#!/bin/python3
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

res = open('../graph.csv', "r")
tmp = res.readline().strip().split(',')
N = int(tmp[0])
M = int(tmp[1])
H = float(tmp[2])
T = float(tmp[3])
m1 = int(tmp[4])
m2 = int(tmp[5])
res.close()
res = np.loadtxt('../graph.csv', delimiter=',', skiprows=1)

tau = T / M
h = H / N
#m1 = math.floor(t1 / tau)
#m2 = math.ceil(t2 / tau)
#print(N, M, H, T, t1, t2, m1, m2)
n = 8
wh = np.linspace(0, H, N + 1, endpoint=True)
#colors = cm.rainbow(np.linspace(0, 1, n + 1))
colors = [
 [5.00000000e-01, 0.00000000e+00, 1.00000000e+00, 1.00000000e+00],
 [2.52941176e-01, 9.25637660e-01, 8.30184031e-01, 1.00000000e+00],
 [2.49019608e-01, 3.84105749e-01, 9.80634770e-01, 1.00000000e+00],
 [1.96078431e-03, 7.09281308e-01, 9.23289106e-01, 1.00000000e+00],
 [1.00000000e+00, 1.22464680e-16, 6.12323400e-17, 1.00000000e+00]
]
#colors = ['0.7', '0.6', '0.5', '0.4', '0']

fig, ax = plt.subplots()
k = 2

left = (m1 + M / 25)
mid = ((m1 + m2) / 2)
right = (m2 - M / 25)

m10 = math.floor(left)
m20 = math.floor(mid)
m30 = math.floor(right)

whf = h * np.arange(N + 1)

ax.plot(whf, res[(m10 - m1) * (N + 1):(m10 - m1 + 1) * (N + 1), 0], color=colors[0], label='v=v_0')
ax.plot(whf, res[(m10 - m1) * (N + 1):(m10 - m1 + 1) * (N + 1), 1], color=colors[1], label=f'v=v_res(t={left*tau:.2})')
ax.plot(whf, res[(m20 - m1) * (N + 1):(m20 - m1 + 1) * (N + 1), 1], color=colors[2], label=f'v=v_res(t={mid*tau:.2})')
ax.plot(whf, res[(m30 - m1) * (N + 1):(m30 - m1 + 1) * (N + 1), 1], color=colors[3], label=f'v=v_res(t={right*tau:.2})')
ax.plot(whf, res[(m10 - m1) * (N + 1):(m10 - m1 + 1) * (N + 1), 2], color=colors[4], linestyle='--', label='v=v_sol')

ax.legend(fontsize=10)
plt.xlabel('x', fontsize=10)
plt.ylabel('x * v', fontsize=10)
#fig.savefig('resgraph.eps', format='eps')
fig.savefig('resgraph.pdf', bbox_inches='tight', format='pdf', transparent=True)
