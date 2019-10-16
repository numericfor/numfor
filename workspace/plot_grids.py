#! /usr/bin/ipython
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt


N = 10
xlin = np.r_[np.linspace(0, 0.5, N, endpoint=False), np.linspace(0.5, 1.4, N)]
xlog = np.geomspace(1.e-5, 10, 20)
plt.ion()
plt.style.use('presentation')

fig = plt.figure(figsize=(10, 7))

plt.plot(xlin, 'o', label='lineally-spaced')
plt.xlabel('index')
plt.xticks([0, 5, 10, 15, 20])
plt.legend()
plt.savefig('figs/ex_linspace.png')


fig = plt.figure(figsize=(10, 7))

plt.plot(xlog, 'o', label='log-spaced')
plt.xlabel('index')
plt.yscale('log')
plt.xticks([0, 5, 10, 15, 20])
plt.legend()
plt.savefig('figs/ex_logspace.png')
