% System Dynamics and Control                             December 27, 2019
% ME 3030-21 WI19
% Prof. M Diaz-Maldonado
%
% Synopsis:
% Obtains the pulse response of an underdamped second-order system by
% applying the superposition property. The weight of the pulse A = 1.
% It's seen that as T -> 0 the response converges to the impulse response.
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

% defines the second-order constants which specify a second-order system:
DR = 0.2;		% damping ratio
Wn = 1.0;		% natural frequency
T = 0.5;		% pulse duration (sec)
P = 1 / T;		% pulse intensity (note the pulse weight is A=1)

% real and imaginary parts of the roots of the characteristic equation
alpha = -DR * Wn;
beta  = Wn * sqrt(1 - DR^2);


% calculates the transient response, y(t), Equation (7.72)
t = linspace(0, 30, 1024);
y = @(t) 1 - exp(alpha * t) .* ( cos(beta*t) - alpha/beta * sin(beta*t) );
Y = P * ( y(t) - y(t - T) .* (t > T) );

yim = beta * (1 + (alpha / beta)^2) * exp(alpha * t) .* sin(beta * t);
plot(t, Y, '-k', 'linewidth', 2);    hold on	% pulse
plot(t, yim, '--', 'linewidth', 2)		% impulse
legend('pulse', 'impulse')

ylim([-0.5, 1.0])
xlabel('time, t, sec')
ylabel('transient response, y(t)')
title('Pulse Response of an Underdamped Second Order Dynamic System')
grid on

print('pulse-response-second-order-system.png','-r600','-dpng')
