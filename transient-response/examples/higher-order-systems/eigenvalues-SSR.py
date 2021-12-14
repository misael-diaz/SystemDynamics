#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
System Dynamics and Control
ME 3030-21 WI19
Prof. M Diaz-Maldonado                                    November 22, 2021

Synopsis:
Shows how to use numpy's linear algebra method ``eig()'' to obtain the
eigenvalues of the system matrix of the SSR of a linear dynamic system.
It's shown that the eigenvalues are identical to the roots of the
characteristic equation.


Copyright (c) 2021 Misael Diaz-Maldonado
This file is released under the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.


References:
[0] CA Kluever, Dynamic Systems: Modeling, Simulation, and Control.
[1] A Gilat, MATLAB: An Introduction with Applications, 6th edition.
[2] numpy.org/doc/stable/reference/generated/numpy.sort_complex.html
"""

from numpy import array
from numpy import roots
from numpy import sort_complex as sort
from numpy.linalg import eig

# gets the roots of the characteristic equation
C = [1, 3, 8, 12]
r = roots(C)
r = sort(r)

# gets the eigenvalues via numpy's linear algebra method `eig()'
A = array([
    [  0,  1,  0],
    [  0,  0,  1],
    [-12, -8, -3]
])

eigs, eigvec = eig(A)
eigs = sort(eigs)

zeroes = r - eigs
print(zeroes)


"""
NOTES:
The negligibly small differences between the roots and the eigenvalues
seen in the `zeroes' array may be attributed to machine precision errors.
"""


"""
Miscellaneous:
real = r[isreal(r)]                         # extracts real root
cplx = r[[not real for real in isreal(r)]]  # extracts complex conjugate

obtains the real and imaginary components of the complex conjugate:
re, im = ( 0.5 * cplx.sum(), abs( 0.5 * (cplx[1] - cplx[0]) ) )
"""
