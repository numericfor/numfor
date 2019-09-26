#! /usr/bin/ipython
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import splrep

pyfname = 'data/pysplrep{:.2f}.dat'
fofname = 'data/testsplrep{:.2f}.dat'
ss = [0, 0.25, 0.5, 1]


def create_data():
  N = 6
  x = np.linspace(0, np.pi, N)
  y = np.sin(x)

  for s in ss:
    tck = splrep(x, y, s=s)
    np.savetxt(
        fname.format(s),
        np.c_[
            tck[0],
            tck[1]],
        fmt="%.14g",
        header=' t           c')


# create_data()
for s in ss:
  tp, cp = np.loadtxt(pyfname.)
