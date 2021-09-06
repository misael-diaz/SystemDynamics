% System Dynamics and Control                             December 27, 2019
% ME 3030-21 WI19
% Prof. M Diaz-Maldonado
%
% Synopsis:
% Reproduces the step response of example 7.8 from Kluever's textbook.
%
%
% Copyright (c) 2021 Misael Diaz-Maldonado
% This file is released under the GNU General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
%
% References:
% [0] CA Kluever, Dynamic Systems: Modeling, Simulation, and Control
% [1] A Gilat, MATLAB: An Introduction with Applications, 6th edition

clear
close all
clc

% mechanical system parameters
J = 0.2;			% moment of inertia
b = 1.6;			% friction coefficient
k = 65.;			% stiffness

% magnitude of step input
A = 2.5;

% standard coefficients of the second-order system
a1 = b / J;	a0 = k / J;	b0 = 1 / J;


% defines the second-order constants which specify a second-order system:

Wn = sqrt(a0);			% natural frequency
DR = a1 / (2 * Wn);		% damping ratio

% real and imaginary parts of the roots of the characteristic equation
alpha = -DR * Wn;		% real
beta  = Wn * sqrt(1 - DR^2);	% imaginary


% calculates the transient response, y(t), Equation (7.72)
t = linspace(0, 2, 512);
y = (  1 - exp(alpha * t) .* ( cos(beta*t) - alpha/beta * sin(beta*t) )  );
y = (b0 * A / Wn^2) * y;
plot(t, y, '-k')
hold on



% specifies the figure labels, the title, and the legend
xlabel('time, t, sec')
ylabel('underdamped step response, y(t)')
title('step response of an underdamped mechanical system')
grid on

% exports figure to PNG format with 600 DPI
print('step-response-mechanical-system.png','-r600','-dpng')
