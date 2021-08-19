% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Course: System Dynamics and Control                       ME 3030-21 WI19
% Author: Prof. M Diaz-Maldonado                          Date: Dec/27/2019
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Synopsis:
% Plots the Bode Diagram for a Second-Order Dynamic System (torsional
% mechanical system). 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
clear;
close all;
clc;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% mechanical system constants:
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% moment of inertia, kg m^2
J = 0.2;
% stiffness constant, N m/rad
k = 65;
% viscous friction, N m s/rad
b = 1.6;


% Constants of the Second Order LTI Ordinary Differential Equation,
% which describes the dynamics of the mechanical system:
%                y'' + a1 y' + a0 y = b0 u(t)
a1 = b/J;
a0 = k/J;
b0 = 1/J;


% natural frequency
omega_n = sqrt(a0);
% damping ratio
DR = a1 / (2*sqrt(a0));
if (DR>0 && DR<1)
    % calculates the resonant frequency if underdamped (0 < DR < 1)
    omega_r = omega_n * sqrt(1-2*DR^2);
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% input u(t):
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% amplitude of input sinusoid, u(t) = U0 sin(wt)
U0 = 1.5; % N m 
% input frequency
w = logspace(-1,4,500);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% frequency response constants, y_ss(t) = A cos(wt) + B sin(wt)
A = -b0 * U0 * (a1*w) ./ ( (a1*w).^2 + (a0-w.^2).^2 ); 
B =  b0 * U0 * (a0 - w.^2) ./ ( (a1*w).^2 + (a0-w.^2).^2 );

% magnitude of amplitude ratio, |G(jw)| = Y0 / U0
% response amplitude, K = Y0
K = sqrt( A.^2 + B.^2 );
% phase angle, < G(jw), radians
phi = -acos(B./K); 
% ----------------------------------------------------------------------- %
% NOTE:
% -acos fun is used to obtain the phase angle, which
% is in the range (0,-pi); other inverse trigonometric
% functions such as asin or atan are constrained to other ranges.
% ----------------------------------------------------------------------- %

% magnitude of amplitude ratio, |G(jw)| = K / U0 = Y0 / U0 
AR = K/U0;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% plots the Bode Diagram
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

subplot(2,1,1);
semilogx(w, 20*log10(AR) );
ylabel('Magnitude of G(jw), dB');
xlabel('Frequency, \omega, rad/s');
title('Bode Diagram for a Second-Order Dynamic System');
grid on;

subplot(2,1,2);
semilogx(w, 180 * phi / pi);
ylabel('Phase of G(jw), deg');
xlabel('Frequency, \omega, rad/s');
grid on;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %