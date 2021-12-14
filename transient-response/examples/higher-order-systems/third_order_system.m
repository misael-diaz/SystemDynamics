% System Dynamics and Control                             December 27, 2019
% ME 3030-21 WI19
% Prof. M Diaz-Maldonado
%
% Synopsis:
% Shows how to use the built-in function roots to obtain the roots of the
% characteristic equation of a dynamic system:
% 	a0 * y''' + a1 * y'' + a2 * y' + a3 * y = b0 * u(t),
%
% where the a's and b's are time-invariant coefficients and u(t) is any
% input function.
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

% coefficients of characteristic equation [a0, a1, a2, a3]
C = [1 2.5 38 18.5];
% obtain the roots of the characteristic equation, stored in vector r
r = roots(C)
