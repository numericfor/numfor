#! /usr/bin/ipython
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import splrep, splev

# pyfname = 'data/pysplrep_s{}.dat'
# fofname = 'data/testsplrep_s{}.dat'
pyfname = 'data/pysplrep_s{}.dat'
fofname = 'data/testsplrep_per_s{}.dat'
pyfdata = 'data/pysplev_s{}.dat'
fofdata = 'data/fosplev_s{}.dat'

ss = [0, 1, 2, 3]
N = 6
x = np.linspace(0, np.pi, N)
y = np.sin(x)
xnew = np.linspace(0, np.pi, 50)

plt.ion()


def create_sin_data():
  for s in ss:
    tck = splrep(x, y, s=s)
    np.savetxt(pyfname.format(s), np.c_[tck[0], tck[1]], fmt="%.14g",
               header=' t           c')
    return tck


def create_sin_spl():
  for s in ss:
    tck = splrep(x, y, s=s)
    ynew = splev(xnew, tck)
    np.savetxt(pyfname.format(s), np.c_[tck[0], tck[1]], fmt="%.14g",
               header=' x           y')


def create_sin_data_per():
  for s in ss:
    tck = splrep(x, y, s=s, per=True)
    np.savetxt(pyfname.format(s), np.c_[tck[0], tck[1]], fmt="%.14g",
               header=' t           c')


# Plot


def plot_sin_t():

  for s in ss:
    figt, axt = plt.subplots(num='t' + str(s))

    tp, cp = np.loadtxt(pyfname.format(s), unpack=True)
    tf, cf = np.loadtxt(fofname.format(s), unpack=True)
    axt.plot(tp, 'o', color='C1')
    axt.plot(tf, 'x', color='C2')

  plt.legend(loc='best')


def plot_sin_c():

  for s in ss:
    figc, axc = plt.subplots(num='c' + str(s))

    tp, cp = np.loadtxt(pyfname.format(s), unpack=True)
    tf, cf = np.loadtxt(fofname.format(s), unpack=True)

    axc.plot(cp, 'o', color='C1', label='py')
    axc.plot(cf, 'x', color='C2', label='fo')
  plt.legend(loc='best')


# create_sin_data_per()


# plt.plot(x, y, 'o')
# plt.plot(xnew, ynew)

# plot_sin_c()
