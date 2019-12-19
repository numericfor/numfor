#! /usr/bin/ipython
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
plt.ion()
x = np.linspace(0, 10, 500)
plt.style.use('presentation')
fig = plt.figure(figsize=(10, 7))
plt.plot(x, np.abs(np.log(x + 1) / (1 + 100 * x**2)))
plt.text(4, 0.03, r"$\frac{\log(x+1)}{1 + 100 x^2}$", fontsize='x-large')
plt.savefig("figs/ex_integr_iadsi.png")
