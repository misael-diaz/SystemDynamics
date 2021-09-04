% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Course: System Dynamics and Control                       ME 3030-21 WI19
% Author: Prof. M Diaz-Maldonado                          Date: Dec/27/2019
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Synopsis:
% Reproduces the response of example 7.8 from Kluever's textbook.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
clear;
close all;
clc;

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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% defines the second-order constants which specify a second-order system:

% natural frequency
Wn = sqrt(a0);
% damping ratio
DR = a1 / (2*Wn);

% real and imaginary parts of the roots of the characteristic equation
alpha = -DR * Wn;
beta  = Wn * sqrt(1-DR^2);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% calculates the transient response, y(t), Equation (7.72)
t = linspace(0,2,500);
y = b0 * A / Wn^2 *...
    (1 - exp(alpha * t) .* (cos(beta*t) - alpha/beta * sin(beta*t)));
plot(t,y,'-k'); hold on;



% specifies the figure labels, the title, and the legend
xlabel('time, t, sec');
ylabel('underdamped unit-step response, y(t)');
title('step response of an underdamped mechanical system');
%title('step response of an underdamped second order system (\omega_n = 1 rad/s)');
grid on;

% exports figure to PNG format with 600 DPI
print('step-response-mechanical-system.png','-r600','-dpng');