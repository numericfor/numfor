#! /usr/bin/ipython
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

plt.ion()

fig = plt.figure(figsize=(10, 7))
X = np.loadtxt('data/rand_gauss_hist.dat')
plt.bar(X[:, 0], X[:, 1], label='rand')
# plt.plot(X[:, 0], X[:, 2], 'o', label='pdf')
plt.legend()
# plt.savefig('figs/ex_gauss_hist.png')
