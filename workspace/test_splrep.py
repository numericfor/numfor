#! /usr/bin/ipython
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import splrep, splprep, splev

# pyfname = 'data/pysplrep_s{}.dat'
# fofname = 'data/testsplrep_s{}.dat'
fofname = 'data/testsplrep_per_s{}.dat'
fofdata = 'data/fosplev_s{}.dat'


plt.ion()


def create_splrep_data(ss):
  "Create data for comparison with numfor splrep & splev"
  pyfname = 'data/pysplrep_s{}.dat'
  pyfdata = 'data/pysplev_s{}.dat'

  N = 6
  Nnew = 59
  x = np.linspace(0, np.pi, N)
  y = np.sin(x)
  xnew = np.linspace(0, np.pi, Nnew)

  for j, s in enumerate(ss):
    tck = splrep(x, y, s=s)
    ynew = splev(xnew, tck)
    head = 't c for s={}'.format(s)
    np.savetxt(pyfname.format(s), np.c_[tck[0], tck[1]], fmt="%.14g", header=head)
    head = 'x   y  for s={}'.format(s)
    np.savetxt(pyfdata.format(s), np.c_[xnew, ynew], fmt="%.14g", header=head)


def create_splprep_data(ss):
  "Create data for comparison with numfor splprep & splev"
  pyfname = 'data/pysplprep_s{}.dat'
  pyfdata = 'data/pysplpev_s{}.dat'

  phi = np.linspace(0, 2. * np.pi, 40)
  r = 0.5 + np.cos(phi)         # polar coords
  x, y = r * np.cos(phi), r * np.sin(phi)    # convert to cartesian

  for s in ss:
    tck, u = splprep([x, y], s=s)
    np.savetxt(pyfname.format(s), np.c_[tck[1][0], tck[1][1]], fmt="%.14g",
               header=' c[0]           c[1]')
    ynew = splev(u, tck)
    head = 'u    x        y ,  for s={}'.format(s)
    np.savetxt(pyfdata.format(s), np.c_[u, ynew[0], ynew[1]], fmt="%.14g", header=head)


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
ss = [0., 0.3]
create_splprep_data(ss)

# plt.plot(x, y, 'o')
# plt.plot(xnew, ynew)

# plot_sin_c()
