% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Course: System Dynamics and Control                       ME 3030-21 WI19
% Author: Prof. M Diaz-Maldonado                          Date: Dec/28/2019
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Synopsis:
% Example 10-1: Transient Response of a DC Motor. One may change the
% system parameters at will to modify the characteristics of the transient
% response, that is the system might exhibit overdamped or underdamped
% dynamics. Code reproduces the results of Kluever's textbook with the
% default values of the system constants.
%
% The damping ratio, the natural frequency, the roots of the characteristic
% equation, the steady-state angular velocity and the DC Gain are computed
% for this second order dynamic system.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
clear;
close all;
clc;

% defines system constants:
% inductance, Henries, H
L = 1.5e-3;
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

% second-order dynamic system constants
a1 = (J * R + b * L) / (J * L);
a0 = (b * R + Kb * Km) / (J * L);

b0 = Km/(J*L);

% natural frequency, rad/s
omega_n = sqrt(a0)
% damping ratio
DR = a1 / (2*sqrt(a0))

% determine the system's characteristic time, tau
if (DR >= 1)
    % if the system is critically damped or overdamped
    
    % roots of characteristic equation obtained by quadratic equation
    r1 = 0.5 * (-a1 + sqrt(a1^2 - 4 * a0))
    r2 = 0.5 * (-a1 - sqrt(a1^2 - 4 * a0))
    % system time for an overdamped system (two distinct negative roots)
    tau = max( -1/r1, -1/r2 );
    
    % settling time, seconds, s
    ts = 4 * tau
    
    t = linspace(0, 5*ts, 1e4);
    y = b0 * A / omega_n^2 * ( 1 + (r2/r1) / ( 1 - r2/r1 ) * exp(r1*t) -...
        1/(1-r2/r1) * exp(r2*t) );
    plot(t, y);
    xlabel('Time, sec');
    ylabel('Transient Response, y(t) = \omega(t), rad/s');
    title('Overdamped Second-Order Dynamic System Response');
    grid on;
    xlim([0 5*ts]);
    
else
    % otherwise, for an underdamped system the characteristic time comes 
    % from the real part of the complex conjugate roots:
    tau = 1/(DR * omega_n);
    
    alpha = -(DR * omega_n);
    beta = omega_n * sqrt(1-DR^2);
    
    % settling time, seconds, s
    ts = 4 * tau
    t = linspace(0, 5*ts, 1e4);
    y = (b0 * A / omega_n^2) * ( 1 - exp(alpha * t) .* ( cos(beta*t) -...
        alpha/beta * sin(beta*t) ) );
    plot(t, y);
    
    xlabel('Time, sec');
    ylabel('Transient Response, y(t) = \omega(t), rad/s');
    title('Underdamped Second-Order Dynamic System Response');
    grid on;
    xlim([0 5*ts]);
end

% settling time, seconds, s
ts = 4 * tau;

% the DC Gain,  (rad/s) / V 
DCG = (b0/a0)

% the steady-state motor speed (derivatives equal to zero), rad/s
w_ss = DCG * A

