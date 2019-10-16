#! /usr/bin/ipython
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt


def plot_lin_log():
  N = 10
  xlin = np.r_[
      np.linspace(
          0, 0.5, N, endpoint=False), np.linspace(
          0.5, 1.4, N)]
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


fig = plt.figure(figsize=(10, 7))
X = np.loadtxt('data/ex_loglinspace.dat')
plt.plot(X[:, 0], 'o', label='step = 1/8')
plt.plot(X[:, 1], 'o', label='step = 1/4')
plt.plot(X[:, 2], 'o', label='step = 1/2')
plt.plot(X[:, 3], 'o', label='step = 1')
plt.legend()
plt.savefig('figs/ex_loglinspace.png')
