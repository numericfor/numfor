#! /usr/bin/ipython
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

plt.ion()

fig = plt.figure(figsize=(10, 7))
# X = np.loadtxt('data/rand_gauss_hist.dat')
X = np.loadtxt('data/rand_exponential_hist.dat')
plt.plot(X[:, 0], X[:, 2], 'C1', lw=4, label='pdf')
plt.plot(X[:, 0], X[:, 1], 'C3', linestyle='--', lw=3, label='rand')
plt.bar(X[:, 0], X[:, 1], width=0.2, alpha=0.7)
plt.legend()
# plt.savefig('figs/ex_gauss_hist.png')
