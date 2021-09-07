% System Dynamics and Control		                  December 27, 2019
% ME 3030-21 WI19
% Prof. M Diaz-Maldonado
%
% Synopsis:
% Obtains the unit-step response of an underdamped second-order system.
% (Reproduces Figure 7.19 of Kluever's textbook.)

clear
close all
clc

% defines the second-order constants which specify a second-order system:
DR = 0.2;			% damping ratio
Wn = 1.0;			% natural frequency

% real and imaginary parts of the roots of the characteristic equation
alpha = -DR * Wn;		% real
beta  = Wn * sqrt(1 - DR^2);	% imaginary

% calculates the transient response, y(t), Equation (7.72)
t = linspace(0, 25, 512);
y = 1 - exp(alpha * t) .* ( cos(beta * t) - alpha / beta * sin(beta * t) );
plot(t, y, '-k')
hold on


% plots the step response for a higher damping ratio on the same figure
DR = 0.4;
alpha = -DR * Wn;
beta  = Wn * sqrt(1 - DR^2);
y = 1 - exp(alpha * t) .* ( cos(beta * t) - alpha / beta * sin(beta * t) );
plot(t, y, '-r')


% specifies the figure labels, the title, and the legend
xlabel('time, t, sec')
ylabel('underdamped unit-step response, y(t)')
title('step response of an underdamped second order system')
legend('damping ratio \zeta = 0.2', 'damping ratio \zeta = 0.4')
grid on

% exports figure to PNG format with 600 DPI
print('unit-step-response-second-order-system.png','-r600','-dpng')
