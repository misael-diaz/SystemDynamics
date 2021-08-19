% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Course: System Dynamics and Control                       ME 3030-21 WI19
% Author: Prof. M Diaz-Maldonado                          Date: Dec/28/2019
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Synopsis:
% Example 10-6-B:
% Closed-loop response of an underdamped PD-controlled mass-damper system.
% The system reaches steady state in less than a second! 
% Code reproduces Figure 10.20 of Kluever's textbook.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
clear;
close all;
clc;

% reference signal, step-input position, m 
A = 0.1;

% plot the reference step input
figure(1);
t = linspace(0, 1.5, 1e4);
% step input, u(t) = x_ref(t) = 0.1 U(t)
u = A * (t>=0);
plot(t, u, '--k'); hold on;

% Proportional Gain, V/m 
Kp = 16;
% Derivative Gain, V s/m
Kd = 4;
    
% second-order dynamic system constants:
%                y'' + a1 y' + a0 y = b0 u(t)
%
a1 = 0.3 + 2 * Kd;
a0 = 2*Kp;
% non-homogeneous constants
b1 = 2 * Kd;
b0 = 2*Kp;

% natural frequency, rad/s
omega_n = sqrt(a0);
% damping ratio
DR = a1 / (2*sqrt(a0));
    
if (DR>0 || DR < 1)
    % plots the response, y(t) for a underdamped system
    
    % real part of the characteristic roots
    alpha = -(DR * omega_n);
    % imaginary part of the characteristic roots
    beta = omega_n * sqrt(1-DR^2);
    
    % step response
    y1 = (b0 * A / a0) * ( 1 - exp(alpha * t) .* ( cos(beta*t) -...
        alpha/beta * sin(beta*t) ) );
    
    % impulse response
    y2 = (b1 * A / a0) * ( beta * (1 + (alpha/beta)^2) ) * ...
        exp(alpha * t) .* sin(beta*t);
    
    % closed-loop response, y(t) = y1(t) + y2(t)
    y = y1+y2;
    
    plot(t, y); hold on;
end
 
% system time
tau = -1/alpha;
ts = 4*tau;

% the overshoot for a pure step input is the following:
% Mos = exp( -pi * DR / sqrt(1-DR^2) );
% note that for this problem there are two inputs: step and impulse


% final figure options
xlabel('Time, sec');
ylabel('Closed-loop Response, y(t) = x(t), m');
title('Response of the PD-Controlled Mass-Damper System');
grid on;

ylim([0, 0.14]);
legend('step-input, x_{ref}(t)', 'K_P = 16 V/m, K_D = 4 V \cdot s/m');