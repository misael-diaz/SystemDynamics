#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
System Dynamics and Control                                 August 26, 2021
Prof. M Diaz-Maldonado


Synopsis:
Obtains the transient response of an underdamped second-order system to a
unit-step input. (Reproduces Figure 7.19 of Kluever's textbook.)


Copyright (c) 2021 Misael Diaz-Maldonado
This file is released under the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.


References:
[0] CA Kluever, Dynamic Systems: Modeling, Simulation, and Control
[1] R Johansson, Numerical Python: Scientific Computing and Data
    Science Applications with NumPy, SciPy, and Matplotlib, 2nd edition
"""

from numpy import linspace
from numpy import cos, sin, exp, sqrt
import matplotlib as mpl
mpl.use("qt5agg")
mpl.rcParams['text.usetex'] = True      # uses TeX for Math Symbols
from matplotlib import pyplot as plt


def step(zeta, omega, t):
    """
    Synopsis:
    Computes the step response y(t) of an underdamped second-order system
    given the damping ratio zeta, the natural frequency omega, and time t.
    """
    # obtains the real and imaginary components of the characteristic roots
    alpha, beta = ( -zeta * omega, omega * sqrt(1 - zeta**2) )

    y = 1 - exp(alpha * t) * ( cos(beta * t) - alpha / beta * 
        sin(beta * t) )
    
    return y


""" defines system parameters """
zetas = [0.2, 0.4]          # damping ratios
omega = 1.0                 # natural frequency, radians / second
t = linspace(0, 25, 500)    # time, seconds


""" plots the response for varying damping ratios on the same figure """
# defines linestyles and labels
colors = ["black", "red"]
linestyles = ["-", "--"]
labels = [r"$\zeta = 0.2$", r"$\zeta = 0.4$"]

plt.close("all")
plt.ion()
fig, ax = plt.subplots()

for n, zeta in enumerate(zetas):
    ax.plot(t, step(zeta, omega, t), color = colors[n],
            linestyle = linestyles[n], label = labels[n])

ax.set_xlabel('time, t, seconds')
ax.set_ylabel('transient response, y(t)')
ax.set_title('Step Response of an Underdamped Second Order Dynamic System')
ax.legend()

fig.savefig("unit-step-response-second-order-dynamic-system.png", dpi=300)
