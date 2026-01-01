#!/bin/python3
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

filename = './graph_partial_solution.csv'
res = open(filename, "r")
tmp = res.readline().strip().split('\t')
H = float(tmp[0])
N = int(tmp[1])
Ns = int(tmp[2])
T = float(tmp[3])
M = int(tmp[4])
res.close()
res = np.loadtxt(filename, delimiter='\t', skiprows=1)


file_analytic = './graph_partial_analytic.csv'
res_analytic = np.loadtxt(file_analytic, delimiter='\t', skiprows=1)

tau = T / M
h = H / N

plt.plot(res[:, 0], res[:, 1], color='blue', label='solution')
plt.plot(res_analytic[:, 0], res_analytic[:, 1], color='red', linestyle='--',  label='analytic')
plt.legend(fontsize=10)
plt.grid()
plt.xlabel(r'$\xi$', fontsize=10)
plt.ylabel(rf'$\Phi(\xi, {T})$', fontsize=10)
plt.savefig('graph_partial.pdf', format='pdf')

plt.clf()

plt.plot(res[:, 0], [x * c for x, c in zip(res[:, 0], res[:, 1])], color='blue', label='solution')
plt.plot(res_analytic[:, 0], [x * c for x, c in zip(res_analytic[:, 0], res_analytic[:, 1])], color='red', linestyle='--',  label='analytic')
plt.legend(fontsize=10)
plt.grid()
plt.xlabel(r'$\xi$', fontsize=10)
plt.ylabel(rf'$\xi\cdot\Phi(\xi, {T})$', fontsize=10)
plt.savefig('graph_partial_m.pdf', format='pdf')
