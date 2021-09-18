#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
System Dynamics and Control                              September 18, 2021
Prof. M. Diaz-Maldonado


Synopsis:
Solves for the transient response of a first-order system when subjected
to a pulse input:

                    y' + k * y = b * u(t),  y(0) = 0,

where k is the rate constant, b is the ``forcing'' constant, and u(t)
is either the unit step or the unit pulse input function of duration T.


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


# initial value, rate and forcing constants, and pulse lapse, respectively
yi, k, b, T = (0.0, 1.0, 1.0, 0.5)
# adjusts the signal intensity so that the pulse has a unit magnitude
P = 1 / T
# defines lambdas for the analytic step, pulse, and impulse-responses
step    = lambda t: (yi - b * P / k) * exp(-k * t) + b * P / k
pulse   = lambda t: step(t) - step(t - T) * (t >= T)
impulse = lambda t: -k * (yi - b / k) * exp(-k * t)
# defines the right-hand side RHS of the ODE: dy/dt = f(t, y) as a lambda
odefun  = lambda t, y: (b * P * (t < T) - k * y)


""" solves the transient response via the 4th-order Runge-Kutta Method """
ti, tf = (0.0, 12.0)            # initial time and final integration times
tspan  = (ti, tf)               # time span
t = linspace(ti, tf, 256)
odesol = solve_ivp(odefun, tspan, [yi], method="RK45")
# unpacks the numerical solution
t_scipy, y_scipy = (odesol.t, odesol.y[0, :])


# plots the response
plt.close("all")
plt.ion()
fig, ax = plt.subplots()

ax.plot(t, pulse(t), color="blue", linewidth=2.0, linestyle="--",
        label="pulse")
ax.plot(t, impulse(t), color="red", linewidth=2.0,
        label="impulse")
ax.plot(t_scipy, y_scipy, color="black", marker="v", linestyle="",
        label="scipy's 4th-order Runge-Kutta method")

ax.grid()
ax.legend()
ax.set_xlabel("time, t")
ax.set_ylabel("dynamic response, y(t)")
ax.set_title("Pulse Response of a First Order Dynamic System")



"""
Comments:
Plots the transient response of the first-order dynamic system to a
unit impulse for comparison. You are encouraged to verify that decreasing
the duration of the pulse, T, yields a transient response closer to that
of a unit impulse.
"""
