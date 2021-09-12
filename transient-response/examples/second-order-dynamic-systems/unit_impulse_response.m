% System Dynamics and Control                       	   January 22, 2020
% ME 3030-21 WI19
% Prof. M Diaz-Maldonado
% Synopsis:
% Plots the response of a second order system to a unit impulse input.

clear
close all
clc

% defines the second-order constants which specify a second-order system:
DR = 0.15;		% damping ratio
Wn = 1.00;		% natural frequency

% real and imaginary parts of the roots of the characteristic equation
alpha = -DR * Wn;
beta  = Wn * sqrt(1 - DR^2);

% calculates the transient response, y(t)
% Note: b0 = 1/Wn^2 and A=1.
t = linspace(0, 30, 1024);
y = beta * (1 + (alpha / beta)^2) * exp(alpha * t) .* sin(beta * t);
plot(t, y, '-k'); hold on

% plots the impulse response for a higher damping ratio on the same figure
DR = 0.3;
alpha = -DR * Wn;
beta  = Wn * sqrt(1 - DR^2);
y = beta * (1 + (alpha / beta)^2) * exp(alpha * t) .* sin(beta * t);
plot(t, y, '-r')

% specifies the figure labels, the title, and the legend
xlabel('time, t, sec')
ylabel('underdamped impulse response, y(t)')
title('impulse response of an underdamped second order system')
grid on;

legend('damping ratio \zeta = 0.15', 'damping ratio \zeta = 0.30')

print('impulse-response-second-order-system.png','-r600','-dpng')
