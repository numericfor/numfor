#! /usr/bin/ipython
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import splrep

N = 6

x = np.linspace(0, np.pi, N)
y = np.sin(x)

tck = splrep(x, y)
