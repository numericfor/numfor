#! /usr/bin/ipython
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

plt.ion()
plt.style.use('paper')
fname = "data/ex_interp_cspline1.dat"

Ndim = 20

fig = plt.figure(figsize=(10, 7))
x0 = np.linspace(1.e-6, 20., Ndim)
y0 = np.sin(x0)
x, y, ys = np.loadtxt(fname, unpack=True)
plt.plot(x, y, '-', label='sin(x)')
plt.plot(x, ys, '--', label='cspline')
plt.plot(x0, y0, 'o', label='data points')
plt.legend(loc='best')

plt.savefig('figs/ex_interp_cspline1.png')
