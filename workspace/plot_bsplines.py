#! /usr/bin/ipython
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

plt.ion()
plt.style.use('presentation')


def plot_splrep1():
  fname = "data/ex_interp_splrep1.dat"

  N = 6

  fig = plt.figure(figsize=(10, 7))
  x0 = np.linspace(0, np.pi, N)
  y0 = np.sin(x0)
  # x, y, ys = np.loadtxt(fname, unpack=True)
  x, y = np.loadtxt(fname, unpack=True)
  plt.plot(x, np.sin(x), '-', label='sin(x)')
  plt.plot(x, y, '--', label='spline')
  plt.plot(x0, y0, 'o', label='data points')
  plt.legend(loc='best')

  plt.savefig('figs/ex_interp_splrep1.png')


def plot_splprep1():
  fname = "data/ex_interp_splprep1.dat"

  Nd = 40

  fig = plt.figure(figsize=(10, 7))

  t = np.linspace(0, 2 * np.pi, Nd)
  r = 0.5 + np.cos(t)
  x0 = r * np.cos(t)
  y0 = r * np.sin(t)

  # x, y, ys = np.loadtxt(fname, unpack=True)
  u, x, y = np.loadtxt(fname, unpack=True)
  plt.plot(x, y, '--', label='spline')
  plt.plot(x0, y0, 'o', label='data points')
  plt.legend(loc='best')

  plt.savefig('figs/ex_interp_splprep1.png')


plot_splrep1()
plot_splprep1()
