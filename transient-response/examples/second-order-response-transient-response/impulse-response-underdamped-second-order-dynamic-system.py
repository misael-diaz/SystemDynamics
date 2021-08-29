#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
System Dynamics and Control                                 August 29, 2021
Prof. M Diaz-Maldonado


Synopsis:
Solves for the impulse response of an underdamped mechanical system
described by the second-order Ordinary Differential Equation ODE:

                m * y'' + b * y' + k * y = f * u(t),

where m is the mass, b is the damper friction, k is the spring stiffness,
f is the forcing constant, and u(t) is the impulse of magnitude P.

The code produces a plot qualitatively similar to that shown in figure 7.22
of Kluever's textbook.


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
from numpy import pi, cos, sin, exp, sqrt, arctan
import matplotlib as mpl
mpl.use("qt5agg")
mpl.rcParams['text.usetex'] = True
from matplotlib import pyplot as plt


def impulse(prms, t):
    """
    Synopsis:
    Computes the impulse-response y(t) of an underdamped second-order
    dynamic system given the system mass, friction and stiffness constants,
    and the forcing constant and magnitude P of the impulse u(t).
    """

    # unpacks mechanical system params (inertia, friction, stiffness, etc.)
    m, b, k, f, P = prms
    # defines the standard coefficients of the second-order ODE
    a0, a1, b0 = (k / m, b / m, f / m)
    # determines the natural frequency and damping ratio, respectively
    omega, zeta = (  sqrt(a0),  a1 / ( 2 * sqrt(a0) )  )
    # complains if the system is not underdamped
    is_underdamped(zeta)
    # obtains the real and imaginary components of the characteristic roots
    alpha, beta = ( -zeta * omega, omega * sqrt(1 - zeta**2) )
    
    # computes the system response y(t) to a unit-impulse input
    y = beta * (1 + (alpha / beta)**2) * exp(alpha * t) * sin(beta * t)
    # rescales to obtain the system response to an impulse of magnitude P
    y *= b0 * P / omega**2

    return y


def perf(prms):
    """
    Synopsis:
    Obtains the performance characteristics of a second-order underdamped
    dynamic system.
    """

    # calculates the underdamped mechanical system characteristics
    a0, a1, b0, omega, zeta = underdamped(prms)
    alpha, beta = ( -zeta * omega, omega * sqrt(1 - zeta**2) )
  
    """ peak time and value """
    tp = arctan(-beta / alpha) / ( omega * sqrt(1 - zeta**2) )
    y_max = impulse(prms, tp)

    y_ss = 0.0                                  # steady-state value
    ts = 4 / (zeta * omega)                     # settling time
    T = 2 * pi / ( omega * sqrt(1 - zeta**2) )  # period of oscillation
    N = 2 * sqrt(1 - zeta**2) / (pi * zeta)     # number of cycles

    return (tp, ts, T, N, y_max, y_ss)


def underdamped(prms):
    """
    Synopsis:
    Returns the basic characteristics of an underdamped mechanical system
    given the mass/moment of inertia and friction and stiffness constants.
    """
    # unpacks mechanical system parameters (inertia, friction, stiffness)
    m, b, k, f, P = prms
    # computes the standard coefficients of the second-order ODE
    a0, a1, b0 = (k / m, b / m, f / m)
    # determines the natural frequency and damping ratio, respectively
    omega, zeta = (  sqrt(a0),  a1 / ( 2 * sqrt(a0) )  )
    # complains if not underdamped
    is_underdamped(zeta)
    return (a0, a1, b0, omega, zeta)


def is_underdamped(zeta):
    # complains if the mechanical system is not underdamped

    errMSG = (f"Not an underdamped mechanical system\n" +
        f"              the damping ratio (zeta = {zeta}) is not " +
        f"less than one\n"
    )

    if zeta >= 1.0:
        raise RuntimeError(errMSG)
    return


""" mechanical system parameters """
# moment of inertia, friction, and stiffness
prms = m, b, k, f, P = (0.2, 1.6,  65.0, 1.0, 2.5)
#prms = m, b, k, f, P = (1.0, 8.0, 325.0, 5.0, 2.5) # equivalent set


""" simulation parameters """
# time t, seconds
t = linspace(0, 2, 512)


""" displays performance characteristics """
tp, ts, Period, Ncycles, y_max, y_ss = perf(prms)
performance = (
    f"\n\n"
    f"peak time and value:             {tp}, {y_max}\n"
    f"settling time:                   {ts}\n"
    f"steady-state value:              {y_ss}\n"
    f"Oscillation Period:              {Period}\n"
    f"Number of Oscillation Cycles:    {Ncycles}\n\n"
)
print(performance)


""" transient response visualization """
plt.close("all")
plt.ion()
fig, ax = plt.subplots()
ax.plot(t, impulse(prms, t), color = "black", linestyle = "-",
        linewidth = 2)
ax.set_xlabel(r'time, $t$, seconds')
ax.set_ylabel(r'transient response, $y(t)$')
ax.set_title(
    'Impulse Response of an Underdamped Second Order Dynamic System'
)
ax.grid()

fig.savefig("impulse-response-underdamped-2nd-order-system.png", dpi=300)



"""
Comments:
The performance equations for the settling time, period of oscillation,
and number of oscillation cycles are the same as those for the step
response. The peak time, however, occurs at a different time which
can be determined analytically.

Obtaining the first-derivative of the impulse response, y'(t), and finding
the time where y'(t) is equal to zero yields the sought expression.
"""
