#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
System Dynamics and Control                              September 12, 2021
Prof. M. Diaz-Maldonado


Synopsis:
Solves for the transient response of a first-order system when subjected
to a sinusoid input:

                    y' + k * y = b * u(t),  y(0) = 0,

where k is the rate constant, b is the ``forcing'' constant, and u(t)
is either the sinusoid input function u(t) = sin(omega * t), where
omega is the oscillation frequency of the input signal.


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
from numpy import pi, exp, sin, cos, linspace
from scipy.integrate import solve_ivp
from scipy.optimize import bisect
import matplotlib as mpl
mpl.use("qt5agg")
from matplotlib import pyplot as plt


# initial value, rate and forcing constants, and oscillation frequency
yi, k, b, omega = (0.0, 0.5, 1.0, 1.0)
""" defines lambdas for the analytic sinusoid response, y(t) """
w = omega
A0, A1 = (-w * b / (w * w + k * k), k * b / (w * w + k * k) )
sinusoid = lambda t: ( A0 * (cos(w * t) - exp(-k * t)) + A1 * sin(w * t) )
# defines the right-hand side RHS of the ODE: dy/dt = f(t, y) as a lambda
odefun = lambda t, y: ( b * sin(omega * t) - k * y )


""" solves the transient response via the 4th-order Runge-Kutta Method """

# finds the ``period'' of the transient response y(t)
p = period = bisect(sinusoid, 1.5 * pi, 2.5 * pi)
ti, tf = tspan = (0.0, p)
t = linspace(ti, tf, 256)
odesol = solve_ivp(odefun, tspan, [yi], method="RK45")
# unpacks the numerical solution
t_scipy, y_scipy = (odesol.t, odesol.y[0, :])


# plots the sinusoid response
plt.close("all")
plt.ion()
fig, ax = plt.subplots()

ax.plot(t, sinusoid(t), color="black", linewidth=2.0,
        label="analytic solution")
ax.plot(t_scipy, y_scipy, color="blue", marker="v", linestyle="",
        label="scipy's 4th-order Runge-Kutta method")

ax.grid()
ax.legend()
ax.set_xlabel("time, t")
ax.set_ylabel("transient response, y(t)")
ax.set_title("Sinusoid Response of a First Order Dynamic System")


"""
Note:
The period of the transient response (the time the system spends executing
a cycle) depends on the characteristics of the system (rate k) and that of
the input signal (frequency omega). Thus, the arguments supplied to the
bisection method to find the period shall need adjusting for other
combinations of those parameters.

If the period of the input signal is large (low frequencies) compared to
the settling time of the system, the transient regime ends before a full
cycle is executed. On the other hand, if the period is small (high
frequencies) the system executes several cycles before the transient
regime comes to an end. (Henceforth the system will continue oscillating
steadily as long the input signal is applied -- frequency response regime.)

In this example the system executes about 1.3 cycles before the transient
regime dies out (using a settling time, ts = 4 / k, for the approximation).
"""
