% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Course: System Dynamics and Control                       ME 3030-21 WI19
% Author: Prof. M Diaz-Maldonado                          Date: Dec/27/2019
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Synopsis:
% Plots the Bode Diagram for a RL-Circuit (First-Order Dynamic System).
% Based on Example 9.1 from Kluever's textbook.
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
% sinusiod input (voltage source): e_in(t) = U0 sin(wt)
% input amplitude, U0, Volts
U0 = 2; 
% RHS constant of the LTI First-Order ODE
b = 1/R;  
% input frequency, rad/s
omega = logspace(-1,4,500); % a logarithmic scale is more useful here
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% characterizes the frequency response (amplitude and phase angle):
% steady-state frequency response solution:
%
%                   y_ss = A cos(wt) + B sin(wt) 
%
% steady-state solution constants (obtained analytically):
A = -tau*omega * b ./ (1 + (tau*omega).^2);
B = b ./ (1 + (tau*omega).^2);
% Amplitude Ratio (Y0/U0 = K/U0)
AR = sqrt(A.^2 + B.^2);
% Phase Angle, phi, for a sinusoid response:    y_ss = K sin(wt + phi).
phi = -atan(tau*omega); 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Plots the Bode Diagram
figure(1);

subplot(2,1,1);
% plots magnitude of the amplitude ratio |G(jw)| with increasing frequency
semilogx(omega, 20*log10(AR) );
ylabel('Amplitude Ratio |G(jw)|, dB');
xlabel('Frequency, \omega, rad/s');
title('Bode Plot of a First-Order Dynamic System (RL-Circuit)');
ylim([-50, 10])

grid on;

subplot(2,1,2);
% plots the phase angle, phi, in degrees with increasing frequency
semilogx(omega, 180 * phi / pi);
ylabel('Phase Angle of G(jw), deg');
xlabel('Frequency, \omega, rad/s');
ylim([-95, 5]);

grid on;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
 