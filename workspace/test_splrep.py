#! /usr/bin/ipython
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import splrep

pyfname = 'data/pysplrep_s{}.dat'
fofname = 'data/testsplrep_s{}.dat'
ss = [0, 0.25, 0.5, 0.75]


def create_sin_data():
  N = 6
  x = np.linspace(0, np.pi, N)
  y = np.sin(x)
  for s in ss:
    tck = splrep(x, y, s=s)
    np.savetxt(pyfname.format(s), np.c_[tck[0], tck[1]], fmt="%.14g",
               header=' t           c')


# create_sin_data()
# Plot
plt.ion()

for s in ss[:]:
  figt, axt = plt.subplots(num='t' + str(s))
  # figc, axc = plt.subplots(num='c' + str(s))

  tp, cp = np.loadtxt(pyfname.format(s), unpack=True)
  tf, cf = np.loadtxt(fofname.format(s), unpack=True)
  axt.plot(tp, 'o', color='C1')
  axt.plot(tf, 'x', color='C2')
  # axc.plot(cp, 'o', color='C1')
  # axc.plot(cf, 'x', color='C2')
