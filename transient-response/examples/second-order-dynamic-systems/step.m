% System Dynamics and Control                             December 27, 2019
% ME 3030-21 WI19
%
% Prof. M Diaz-Maldonado
% Synopsis:
% Reproduces the step response of example 7.8 from Kluever's textbook.
%
clear
close all
clc

% mechanical system parameters
% moment of inertia
J = 0.2;
% friction coefficient
b = 1.6;
% stiffness
k = 65;

% magnitude of step input
A = 2.5;

% standard coefficients of the second-order system
a1 = b/J; a0 = k/J; b0 = 1/J;


% defines the second-order constants which specify a second-order system:

% natural frequency
Wn = sqrt(a0);
% damping ratio
DR = a1 / (2*Wn);

% real and imaginary parts of the roots of the characteristic equation
alpha = -DR * Wn;
beta  = Wn * sqrt(1-DR^2);


% calculates the transient response, y(t), Equation (7.72)
t = linspace(0,2,500);
y = (1 - exp(alpha * t) .* (cos(beta*t) - alpha/beta * sin(beta*t)));
y = (b0 * A / Wn^2) * y;
plot(t, y, '-k'); hold on;



% specifies the figure labels, the title, and the legend
xlabel('time, t, sec');
ylabel('underdamped step response, y(t)');
title('step response of an underdamped mechanical system');
grid on;

% exports figure to PNG format with 600 DPI
print('step-response-mechanical-system.png','-r600','-dpng');
