#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
System Dynamics and Control                                 August 29, 2021
Prof. M Diaz-Maldonado


Synopsis:
Solves the Initial Value Problem IVP that models the underdamped response,
y(t), of an underdamped third-order dynamic system:

y''' + a2 * y'' + a1 * y1 + a0 * y = b0 * u(t), y'' = y' = y = 0 at t = 0,

where a0, a1, a2, and b0 are the ``standard'' coefficients of the Ordinary
Differential Equation ODE, and u(t) is the step response of magnitude h.


Copyright (c) 2021 Misael Diaz-Maldonado
This file is released under the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.


References:
[0] CA Kluever, Dynamic Systems: Modeling, Simulation, and Control
[1] R Johansson, Numerical Python: Scientific Computing and Data
    Science Applications with NumPy, SciPy, and Matplotlib, 2nd edition
[2] (linsolve: docs.sympy.org/latest/modules/solvers/solveset.html)
[3] (sequence unpacking: stackoverflow.com/questions/38418205/
     get-a-value-from-solution-set-returned-as-finiteset-by-sympy)
"""


import sympy
from sympy import diff as D
from numpy import pi, sqrt
from numpy import empty, linspace
from scipy.optimize import bisect
from scipy.integrate import solve_ivp
import matplotlib as mpl
mpl.use("qt5agg")
from matplotlib import pyplot as plt


def fstep(R, Alpha, Beta, B0, H):

    # time
    t = sympy.Symbol('t')
    # characteristic roots for the underdamped dynamic system
    r = sympy.Symbol('r')                           # real
    alpha, beta = sympy.symbols( ('alpha', 'beta') )# complex conjugates


    # expresses the standard coefficients of the ODE in terms of the roots
    a0 = -r * (alpha**2 + beta**2)
    a1 = 2 * r * alpha  + (alpha**2 + beta**2)
    a2 = -(r + 2 * alpha)
    b0 = sympy.Symbol('b0')


    # defines the homogeneous solution, y
    y1 = sympy.exp(r * t)
    y2 = sympy.exp(alpha * t) * sympy.sin(beta * t)
    y3 = sympy.exp(alpha * t) * sympy.cos(beta * t)

    y = y1 + y2 + y3

    # defines the particular solution, yp
    h = sympy.Symbol('h')
    u = h                   # u(t), step
    yp = b0 * h / a0




    """ solves the Initial Value Problem IVP """
    # matrix
    Mat = sympy.Matrix([
        [y1.subs({t: 0}), y2.subs({t: 0}), y3.subs({t: 0})],
        [D(y1,t).subs({t: 0}), D(y2,t).subs({t: 0}), D(y3,t).subs({t: 0})],
        [D(y1,t,t).subs({t: 0}), D(y2,t,t).subs({t: 0}),
            D(y3,t,t).subs({t: 0})],
    ])
    # coefficient vector
    vec = sympy.Matrix( [-yp.subs({t: 0}), 0, 0] )

    # obtains the undetermined coefficient of the general solution
    A, B, C = sympy.symbols( ('A', 'B', 'C') )
    """Note: trailing comma enables sequence unpacking of the Finite Set"""
    sol, = sympy.linsolve( (Mat, vec), A, B, C )

    # defines the general solution (transient response `y(t)')
    A, B, C = sol
    y = A * y1 + B * y2 + C * y3 + yp

    step = sympy.lambdify(
        t, y.subs({r: R, alpha: Alpha, beta: Beta, b0: B0, h: H}), 'numpy'
    )

    impulse = sympy.lambdify(
        t, D(y, t).subs({r: R, alpha: Alpha, beta: Beta, b0: B0, h: H}), 'numpy'
    )

    return (step, impulse)


def perf(r, alpha, beta, step, impulse):
    """
    Synopsis:
    Obtains the performance characteristics of the dynamic system.
    """

    omega = sqrt(alpha**2 + beta**2)                    # natural frequency
    zeta  = -alpha / omega                              # damping ratio
    
    ts = max(-4 / r,  -4 / alpha)                       # settling time
    P = Period = 2 * pi / ( omega * sqrt(1 - zeta**2) ) # period
    N = Ncycles = ts / Period                           # number of cycles

    # Note: See comments at the end of the source regarding the peak time
    tp = bisect(impulse, 0.15 * P, 0.85 * P)            # peak time
    y_max = step(tp)                                    # peak value
    
    A0 = ( -R * (Alpha**2 + Beta**2) )                  # A0 = a0
    y_ss = B0 * H / A0                                  # steady-state


    """ displays the performance characteristics of the dynamic system """
    performance = (
        f"\n\n"
        f"Damping Ratio:                   {zeta}\n"
        f"Natural Frequency:               {omega}\n"
        f"peak time and value:             {tp}, {y_max}\n"
        f"Settling time:                   {ts}\n"
        f"steady-state value:              {y_ss}\n" 
        f"Maximum Overshoot:               {y_max / y_ss - 1}\n"
        f"Oscillation Period:              {Period}\n"
        f"Number of Oscillation Cycles:    {Ncycles}\n\n"
    )

    print(performance)

    return (omega, zeta, ts, tp, P, N, y_max)




""" plots the transient response """
time = linspace(0, 25, 256)

R, Alpha, Beta, B0, H = (-1, -0.25, 1, 1, 1)

step, impulse = fstep(R, Alpha, Beta, B0, H)


def odefun(t, y):
    """
    Synopsis:
    Defines the third-order system as an equivalent system of first-order
    ODEs to obtain the transient response y(t) numerically.
    """

    # defines the ``standard'' coefficients of the third-order ODE
    A0 = -R * (Alpha**2 + Beta**2)
    A1 = 2 * R * Alpha  + (Alpha**2 + Beta**2)
    A2 = -(R + 2 * Alpha)

    f = empty(3)

    f[0] = y[1]
    f[1] = y[2]
    f[2] = -A2 * y[2] - A1 * y[1] - A0 * y[0] + B0 * H

    return f


# solves the third-order ODE numerically for verification
tspan = [0, 25]
""" initial values: y(0) = 0, y'(0) = 0, y''(0) = 0 """
yi = [0, 0, 0]
# applies the fourth-order Runge-Kutta method
odesol = solve_ivp(odefun, tspan, yi, method="RK45")
t, y = odesol.t, odesol.y
y_step = y[0, :]    # step response
y_impulse = y[1, :] # impulse response


plt.close("all")
plt.ion()

fig, ax = plt.subplots()
ax.plot(time, step(time), linewidth = 2, color = "black",
        label = "analytic")
ax.plot(t, y_step, linestyle = "", marker = "d", color = "orange",
        label = "numeric")
ax.set_xlabel("time, t")
ax.set_ylabel("transient response, y(t)")
ax.set_title("Step Response of an Underdamped Third-Order Dynamic System")
ax.legend()
ax.grid()


# displays performance characteristics on the console
perf(R, Alpha, Beta, step, impulse)



"""
Comments on the Determination of the Peak Time:
If the dynamics of the third-order system is dominated by the complex
roots (rather than the real root) the peak time can be observed in the
first oscillation cycle as for second-order systems. Only then one can
solve for the peak time by applying a numerical technique. Note that
solving for the peak time entails solving the nonlinear equation (in time)
y'(t) = 0. Here we applied the bisection method, though the bracketing
interval may need adjusting for other third-order underdamped systems
(different set of parameters). You have the option of applying other
numerical techniques such as Newton-Raphson.

If the real root dominates, the system will exhibit weak oscillations
and its transient response will resemble that of a first-order system.
And its maximum value of y(t) is the steady-state value. This means
that the bisection method will complain that there's no solution in
the given interval. You may want to comment out that part of the code.
"""
