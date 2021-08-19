% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Course: System Dynamics and Control                       ME 3030-21 WI19
% Author: Prof. M Diaz-Maldonado                          Date: Dec/27/2019
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Synopsis:
% Plots the Frequency Response of a Second-Order Dynamic System (torsional
% mechanical system). Reproduces figure 9.11 of Kluever's textbook, note
% that the transient dynamics are disregarded in our calculations.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
clear;
close all;
clc;

% physical constants
% moment of inertia, kg m^2
J = 0.2;
% stiffness constant, N m/rad
k = 65;
% viscous friction, N m s/rad
b = 1.6;

% amplitude of input sinusoid, u(t) = U0 sin(wt)
U0 = 1.5; % N m 
% frequency of input sinusoid, rad/sec
w=18; % w=\omega

% Constants of the Second Order LTI Ordinary Differential Equation,
% which describes the dynamics of the mechanical system:
%                y'' + a1 y' + a0 y = b0 u(t)
a1 = b/J;
a0 = k/J;
b0 = 1/J;

% frequency response constants, y_ss(t) = A cos(wt) + B sin(wt)
A = -b0 * U0 * (a1*w) / ( (a1*w)^2 + (a0-w^2)^2 ); 
B =  b0 * U0 * (a0 - w^2) / ( (a1*w)^2 + (a0-w^2)^2 );

% response amplitude, K = U0 | G(jw) |, determined analytically
K = sqrt( A^2 + B^2 );
% phase angle, < G(jw), radians (determined from trigonometric identity)
phi = -acos(B/K);

% time vector
t = linspace(0, 2, 1e3);
% calculates input sinusoid, u(t)
u = U0 * sin(w*t);
% calculates the frequency response, y_ss(t):
%       y_ss(t) = A * cos(w*t) + B * sin(w*t) = K * sin(w*t + phi) 
y_ss = K * sin(w*t+phi);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% plots the input and the frequency response of the same figure, 
% two y-axes are used for clarity since the amplitudes differ greatly.
colororder({'k', 'r'});

yyaxis left; 
plot(t, u); 
% y-axis 1: input torque T_in(t)
ylabel('Input Torque, T_{in}(t), N \cdot m')
ylim([-2 2]);

yyaxis right;
plot(t, y_ss,'--'); 
grid on;
xlabel('Time, s');
% y-axis 2: angle, \theta(t)
ylabel('Angular Position, \theta(t), rad');
title('Frequency Response of a Second-Order Dynamic System');