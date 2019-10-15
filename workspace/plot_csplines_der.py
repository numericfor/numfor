#! /usr/bin/ipython
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

plt.ion()
plt.style.use('paper')
fname = "data/ex_interp_cspline2.dat"

Ndim = 20

fig = plt.figure(figsize=(10, 7))

x, y, ys = np.loadtxt(fname, unpack=True)
plt.plot(x, y, '-', label='cos(x)')
plt.plot(x, ys, '--', label='cspline')
plt.legend(loc='best')

plt.savefig('figs/ex_interp_cspline2.png')
