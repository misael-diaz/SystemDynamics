#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
System Dynamics and Control                               December 14, 2021
ME 3030-21 WI19
Prof. M Diaz-Maldonado

Synopsis:
Shows how to use the function `roots' to obtain the roots of the
characteristic equation of a third-order dynamic system:

        a0 * y''' + a1 * y'' + a2 * y' + a3 * y = b0 * u(t),

where the a's and b's are time-invariant coefficients and u(t) is any input
function.


Copyright (c) 2021 Misael Diaz-Maldonado
This file is released under the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.


References:
[0] CA Kluever, Dynamic Systems: Modeling, Simulation, and Control.
[1] A Gilat, MATLAB: An Introduction with Applications, 6th edition.
[2] R Johansson, Numerical Python: Scientific Computing and Data Science
    Applications with NumPy, SciPy, and Matplotlib, 2nd edition.
"""

from numpy import roots

# time-invariant coefficients [a0,a1,a2,a3] of the characteristic equation
C = [1, 2.5, 38, 18.5]
# obtains the roots of the characteristic equation, stored in vector r
r = roots(C)
