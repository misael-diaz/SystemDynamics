% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Course: System Dynamics and Control                       ME 3030-21 WI19
% Author: Prof. M Diaz-Maldonado                          Date: Dec/28/2019
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Synopsis:
% Example 10-3: Transient Response of a DC Motor (simplified model). 
% Reproduces the results of Example 10-3 of Kluever's textbook and plots
% the transient response (motor angular velocity).
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
clear;
close all;
clc;

% defines system constants:
% inductance, Henries, H
L = 0;
% resistance, Ohms, \Omega
R = 0.5;
% motor-torque constant, N m / A
Km = 0.05;
% back EMF constant, V s/rad
Kb = 0.05;
%  motor's moment of inertia, kg m^2
J = 2.5e-4; % overdamped
%J = 4.e-6; % underdamped
% bearing friction, N m s/rad
b = 1e-4;
% source voltage step-input, Volts, V
A = 8;

% first dynamic system constants
tau = (R * J) / (b * R + Kb * Km);
% non-homogeneous coefficient
b0 = Km/(b * R + Kb * Km);

% settling time, seconds, s
ts = 4 * tau

% the DC Gain,  (rad/s) / V 
DCG = b0

% the steady-state motor speed (derivatives equal to zero), rad/s
w_ss = DCG * A

% plot the first-order response
t = linspace(0, 5*tau, 1e4);
y = b0 * A * ( 1 - exp(-t/tau) );
plot(t, y, '-k');
grid on;
xlabel('Time, sec');
ylabel('Transient Response, y(t) = \omega(t), rad/s');
title('First Order Transient Response of the DC Motor');

xlim([0 5*tau]);
