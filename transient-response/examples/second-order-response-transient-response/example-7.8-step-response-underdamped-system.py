#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
System Dynamics and Control                                 August 27, 2021
Prof. M Diaz-Maldonado


Synopsis:
Solves for the step response of an underdamped mechanical system.
Produces a plot similar to that shown in example 7.8 of Kluever's textbook.


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
from numpy import pi, cos, sin, exp, sqrt
import matplotlib as mpl
mpl.use("qt5agg")
mpl.rcParams['text.usetex'] = True
from matplotlib import pyplot as plt


def step(prms, H, t):
    """
    Synopsis:
    Computes the step-response y(t) of an underdamped second-order dynamic
    system given the system mass, friction and stiffness constants, and the
    magnitude H of the step-input u(t).
    """

    # defines the standard coefficients of the second-order ODE
    a0, a1, b0 = (k / m, b / m, 1 / m)
    # determines the natural frequency and damping ratio, respectively
    omega, zeta = (  sqrt(a0),  a1 / ( 2 * sqrt(a0) )  )
    # complains if the system is not underdamped
    is_underdamped(zeta)
    # obtains the real and imaginary components of the characteristic roots
    alpha, beta = ( -zeta * omega, omega * sqrt(1 - zeta**2) )
    
    # computes the system response y(t) to a unit-step input
    y = 1 - exp(alpha * t) * ( cos(beta * t) - alpha / beta *
        sin(beta * t) )
    # rescales to obtain the system response to a step of magnitude H
    y *= b0 * H / omega**2

    return y


def perf(prms, H):
    """
    Synopsis:
    Obtains the performance characteristics of a second-order underdamped
    dynamic system.
    """

    # calculates the underdamped mechanical system characteristics
    a0, a1, b0, omega, zeta = underdamped(prms)
  
    tp = pi / ( omega * sqrt(1 - zeta**2) )     # peak time
    y_max = step(prms, H, tp)                   # peak value
    y_ss = b0 * H / omega**2                    # steady-state value
    ts = 4 / (zeta * omega)                     # settling time
    P = 2 * pi / ( omega * sqrt(1 - zeta**2) )  # period of oscillation
    N = 2 * sqrt(1 - zeta**2) / (pi * zeta)     # number of cycles
    Mos = exp(-pi * zeta / sqrt(1 - zeta**2) )  # maximum overshoot

    return (tp, ts, P, N, Mos, y_max, y_ss)


def underdamped(prms):
    """
    Synopsis:
    Returns the basic characteristics of an underdamped mechanical system
    given the mass/moment of inertia and friction and stiffness constants.
    """
    # unpacks mechanical system parameters (inertia, friction, stiffness)
    m, b, k = prms
    # computes the standard coefficients of the second-order ODE
    a0, a1, b0 = (k / m, b / m, 1 / m)
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
prms = m, b, k = (0.2, 1.6, 65.0)


""" simulation parameters """
# time t, seconds
t = linspace(0, 2, 256)
# step-input magnitude, u(t) = H for t > 0
H = 2.5


""" displays performance characteristics """
tp, ts, Period, Ncycles, Mos, y_max, y_ss = perf(prms, H)
performance = (
    f"\n\n"
    f"peak time and value:             {tp}, {y_max}\n"
    f"settling time:                   {ts}\n"
    f"steady-state value:              {y_ss}\n"
    f"Maximum Overshoot:               {Mos}\n"
    f"Oscillation Period:              {Period}\n"
    f"Number of Oscillation Cycles:    {Ncycles}\n\n"
)
print(performance)


""" transient response visualization """
plt.close("all")
plt.ion()
fig, ax = plt.subplots()
ax.plot(t, step(prms, H, t), color = "black", linestyle = "-",
        linewidth = 2)
ax.set_xlabel(r'time, $t$, seconds')
ax.set_ylabel(r'transient response, $y(t)$')
ax.set_title('Step Response of an Underdamped Second Order Dynamic System')
ax.grid()

fig.savefig("figure-7.21.png", dpi=300)



"""
Comments:
The ordinary differential equation that describes the dynamics of the
mechanical system is:
                    y'' + a1 * y' + a0 * y = b0 * u(t),
where a1, a0, and b0 are the standard coefficients, u(t) is the input, and
y = y(t) is the system response. The constants are defined for a mechanical
system in the step() and underdamped() functions.
"""
