% System Dynamics and Control
% ME 3030-21 WI19
% Prof. M Diaz-Maldonado                                  December 27, 2019
%
% Synopsis:
% Shows how to use the eig() built-in function to obtain the eigenvalues of
% the system matrix of the SSR of a linear dynamic system. It's shown that
% the eigenvalues are identical to the roots of the characteristic
% equation.
%
%
% Copyright (c) 2021 Misael Diaz-Maldonado
% This file is released under the GNU General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
%
% References:
% [0] CA Kluever, Dynamic Systems: Modeling, Simulation, and Control.
% [1] A Gilat, MATLAB: An Introduction with Applications, 6th edition.

clear
close all
clc

% gets the roots of the characteristic equation
C=[1 3 8 12];
r=roots(C);
% roots are sorted for ease of comparison
r=sort(r)

% gets the system matrix eigenvalues via the built-in function `eig()'
A = [ 0  1  0
      0  0  1
    -12 -8 -3];
eigs=eig(A);
eigs=sort(eigs)

zeroes = r - eigs

% NOTE:
% The negligibly small differences between the roots and the eigenvalues
% seen in the `zeroes' array may be attributed to machine precision errors.
