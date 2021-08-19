% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Course: System Dynamics and Control                       ME 3030-21 WI19
% Author: Prof. M Diaz-Maldonado                          Date: Dec/28/2019
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Synopsis:
% Example 10-5-B: 
% Reproduces Figure 10.17 of Kluever's textbook. Closed-loop response of a
% DC Motor with PI Controller. 
%
% Shows how the settling time of the PI-Controlled DC Motor is affected by
% increasing the integral gain, Ki while keeping constant the proportional
% gain, Kp. The PI-Controlled DC Motor exhibits a transient response of a
% second-order dynamic system with the natural frequency of oscillation
% being determined solely by the integral gain. Increasing the proportional
% gain tends to enhance the damping ratio. Thus, the performance of
% the DC Motor can be tuned by adjusting Kp and Ki. 
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
% reference input (anglular velocity u(t) = A U(t)), rad/s
A = 50;

% proportional gain
Kp = [0.2, 0.5, 1.0];
% integral gain
Ki = [4, 10, 20];

figure(1);
% plot the input
t = linspace(0, 0.08, 1e3);
u = A * (t>=0);
plot(t, u, '--k'); hold on;

for n=1:length(Kp)
    % plot the closed-loop response (second-order overdamped system)
    
    % second-order dynamic system constants
    a1 = ( b * R + Km * (Kb + Kp(n)) ) / (J * R);
    a0 = (Km * Ki(n)) / (J * R);
    
    % non-homogeneous coefficients
    b1 = (Km * Kp(n)) / (R * J);
    b0 = (Km * Ki(n)) / (R * J);
    
    % natural frequency, rad/s
    omega_n = sqrt(a0)
    % damping ratio
    DR = a1 / (2*sqrt(a0))
    if (DR < 1)
        % halt execution if the system is underdamped
        fprintf('underdamped system\n');
        return
    end
    
    % roots of characteristic equation obtained by quadratic equation
    r1 = 0.5 * (-a1 + sqrt(a1^2 - 4 * a0))
    r2 = 0.5 * (-a1 - sqrt(a1^2 - 4 * a0))
    % system time for an overdamped system (two distinct negative roots)
    tau = max( -1/r1, -1/r2 );
    
    % settling time, seconds, s
    ts = 4 * tau;
    
    
    % y(t) = y1(t) + y2(t) by the superposition principle
    y1 = b0 * A / omega_n^2 * ( 1 + ...
        (r2/r1) / ( 1 - r2/r1 ) * exp(r1*t) - 1/(1-r2/r1) * exp(r2*t) );
    
    y2 = b1 * A / omega_n^2 * ( r2 / ( 1 - r2/r1 ) * exp(r1*t) -...
        r2/(1-r2/r1) * exp(r2*t) );
    
    % transient response
    y = y1 + y2;
    plot(t, y);
end

xlabel('Time, s');
ylabel('DC Motor speed, rad/s');
title('Transient Response of PI-Controlled DC Motor');

xlim([0 0.08]);
ylim([0 60])
grid on;

legend('u(t) = \omega_{ref}', 'Kp = 0.2 V \cdot s/rad, Ki = 4 V/rad',...
    'Kp = 0.5 V \cdot s/rad, Ki = 10 V/rad',...
    'Kp = 1.0 V \cdot s/rad, Ki = 20 V/rad', 'location', 'SE');