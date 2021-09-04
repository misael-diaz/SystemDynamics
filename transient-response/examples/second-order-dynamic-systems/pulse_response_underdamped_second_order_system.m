% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Course: System Dynamics and Control                       ME 3030-21 WI19
% Author: Prof. M Diaz-Maldonado                          Date: Dec/27/2019
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Synopsis:
% Obtains the pulse response of an underdamped second-order system by
% applying the superposition property. The weight of the pulse A = 1.
% It's seen that as T -> 0 the response converges to the impulse response.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
clear;
close all;
clc;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% defines the second-order constants which specify a second-order system:
% damping ratio
DR = 0.2;
% natural frequency
Wn = 1;
% pulse duration (sec)
T = 0.5;
% pulse intensity (note the pulse weight is A=1)
P = 1/T;

% real and imaginary parts of the roots of the characteristic equation
alpha = -DR * Wn;
beta  = Wn * sqrt(1-DR^2);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% calculates the transient response, y(t), Equation (7.72)
t = linspace(0,30,1e3);
y = @(t) 1 - exp(alpha * t) .* (cos(beta*t) - alpha/beta * sin(beta*t));
Y = P * ( y(t) - y(t-T) .* (t>T) );

yim = beta * (1+(alpha/beta)^2) * exp(alpha * t) .* sin(beta*t);
plot(t,Y,'-k'); hold on;
plot(t,yim,'--');
legend('pulse response', 'impulse response');

% specifies the figure labels, the title, and the legend
ylim([-0.5, 1.0])
xlabel('time, t, sec');
ylabel('transient response, y(t)');
title('pulse response of an underdamped second order system');
%title('step response of an underdamped second order system (\omega_n = 1 rad/s)');
grid on;

% exports figure to PNG format with 600 DPI
print('pulse-response-second-order-system.png','-r600','-dpng');