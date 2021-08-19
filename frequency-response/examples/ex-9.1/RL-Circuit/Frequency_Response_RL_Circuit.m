% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Course: System Dynamics and Control                       ME 3030-21 WI19
% Author: Prof. M Diaz-Maldonado                          Date: Dec/26/2019
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Synopsis:
% Plots the frequency response of a RL-Circuit (first-order dynamic
% system). Reproduces figure 9.6 of Kluever's textbook.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
clear;
close all;
clc;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% system parameters:
% resistance, Omhs
R = 1.2;
% inductance, Henries 
L = 0.02;
% system time constant (first-order dynamic system), seconds
tau = L/R;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% sinusiod input (voltage source)
% e_in(t) = U0 sin(wt)
% input amplitude, V
U0 = 2; 
b = 1/R; %(input current amplitude) 
% input frequency, rad/s
omega = 50;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% characterizes the frequency response (amplitude and phase angle):
% steady-state frequency response solution:
%
%                   y_ss(t) = A cos(wt) + B sin(wt)
%
% steady-state solution constants (obtained analytically):
A = -tau*omega * (b * U0) / (1 + (tau*omega)^2);
B = (b * U0) / (1 + (tau*omega)^2);
% Amplitude, K, expression obtained from the trigonometric identity:
%
%         A cos(theta) + B sin(theta) = K sin(theta + phi),
% 
% where K = sqrt(A^2 + B^2), and phi = atan(A/B);
%
% NOTE:
% Identity holds for phase angles: -pi/2 < phi < pi/2.
%
% amplitude of frequency response
K = sqrt(A^2 + B^2);
% Phase Angle, phi, for a sinusoid response:    y_ss(t) = K sin(wt + phi)
phi = atan(A/B);  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% plot the sinusoid input and the frequency response on the same figure
t = linspace(0,0.5,500);
f = U0 * sin(omega*t);
y_ss = K * sin(omega*t + phi);
plot(t,f,t,y_ss); hold on;

legend('input voltage, u(t)=e_{in}(t)', 'output current, y(t)=I_{ss}(t)');
xlabel('Time, s');
ylim([-3, 3]);
grid on;
title('Frequency Response of a First-Order Dynamic System')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
 