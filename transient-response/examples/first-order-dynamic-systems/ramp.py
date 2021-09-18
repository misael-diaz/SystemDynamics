#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
System Dynamics and Control                              September 18, 2021
Prof. M. Diaz-Maldonado


Synopsis:
Solves for the transient response of a first-order system when subjected
to a ramp input:

                    y' + k * y = b * u(t),

where k is the rate constant, b is the ``forcing'' constant, and u(t)
is either the unit-ramp input function.


Copyright (c) 2021 Misael Diaz-Maldonado
This file is released under the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.


References:
[0] A Gilat and V Subramanian, Numerical Methods for Engineers and
    Scientists: An Introduction with Applications using MATLAB
[1] R Johansson, Numerical Python: Scientific Computing and Data
    Science Applications with NumPy, SciPy, and Matplotlib, 2nd edition
[2] CA Kluever, Dynamic Systems: Modeling, Simulation, and Control
"""


import numpy as np
from numpy import exp, linspace
from scipy.integrate import solve_ivp
import matplotlib as mpl
mpl.use("qt5agg")
from matplotlib import pyplot as plt


# initial value, and rate and forcing constants, respectively
yi, k, b = (0.0, 1.0, 1.0)
# defines lambda for the analytic ramp response
ramp   = lambda t: ( b / (k * k) * (exp(-k * t) - 1.0) + b / k * t )
# defines the right-hand side RHS of the ODE: dy/dt = f(t, y) as a lambda
odefun = lambda t, y: (b * t - k * y)


""" solves the transient response via the 4th-order Runge-Kutta Method """
ti, tf = (0.0, 6.0)             # initial time and final integration times
tspan  = (ti, tf)               # time span
t = linspace(ti, tf, 256)
odesol = solve_ivp(odefun, tspan, [yi], method="RK45")
# unpacks the numerical solution
t_scipy, y_scipy = (odesol.t, odesol.y[0, :])


# plots the response
plt.close("all")
plt.ion()
fig, ax = plt.subplots()

ax.plot(t, ramp(t), color="black", linewidth=2.0,
        label="analytic solution")
ax.plot(t_scipy, y_scipy, color="blue", marker="v", linestyle="",
        label="scipy's 4th-order Runge-Kutta method")

ax.grid()
ax.legend()
ax.set_xlabel("time, t")
ax.set_ylabel("transient response, y(t)")
ax.set_title("Ramp Response of a First Order Dynamic System")
