% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Course: System Dynamics and Control                       ME 3030-21 WI19
% Author: Prof. M Diaz-Maldonado                          Date: Dec/28/2019
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Synopsis:
% Example 10-5-A: 
% The influence of the gain of the proportional controller on the transient
% response of the DC motor is shown through plots of the transient
% response with increasing P-controller gain, Kp. It's found that
% increasing Kp improves the steady-state tracking and improves the
% performance of the DC Motor by decreasing the settling time; that is, 
% the DC Motor reaches steady-state sooner. This is confirmed by subsequent
% plots of the proportional gain with the settling time and DC Gain. The
% DC Gain approaches unity, which is characteristic of dynamic systems with
% perfect steady-state tracking with respect to the reference signal
% (step-input in this example).
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


figure(1);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% plots transient response of P-Controlled DC-Motor for different values of
% the proportional controller gain Kp.

% time, sec (vector)
t = linspace(0, 0.08, 1e3);
% input
u = A * (t>=0);

% proportional controller gain
Kp = [0.2, 0.5, 1];
% plots the reference angular velocity (input)
plot(t, u, '--k'); hold on;

for i=1:length(Kp)

    % first-order dynamic system constants

    % system time constant
    tau = R * J / ( b * R + Km * (Kb + Kp(i)) );
    % DC Gain, K
    K = Kp(i) * Km / ( b * R + Km * (Kb + Kp(i)) ); b0 = K;

    % transient response, y(t), angular velocity of DC Motor
    y = A * b0  * ( 1 - exp(-t/tau) );
    plot(t, y);

end

% final figure options
ylim([0, 60]);
xlim([0, 0.08]);
grid on;

xlabel('Time, sec');
ylabel('Transient Response, y(t)=\omega(t)');
title('Proportional Gain Effect on the DC Motor Dynamics');
legend('\omega_{ref} = 50 rad/s', 'K_p = 0.2 V \cdot s/rad',...
    'K_p = 0.5 V \cdot s/rad', 'K_p = 1.0 V \cdot s/rad', 'location', 'SE');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %


figure(2);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% plots the settling time with increasing proportional gain, Kp

Kp = linspace(0, 2);
% system time constant
tau = R * J ./ ( b * R + Km * (Kb + Kp) );
% settling time, ts
ts = 4 * tau;

plot(Kp, ts, '-r');

xlabel('Gain of Proportional Controller, Kp');
ylabel('Settling time, ts');
title('P-Controller Gain influence over the settling time');

grid on;

figure(3);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% plots the DC Gain with increasing proportional gain, Kp

% DC Gain, K
K = Kp * Km ./ ( b * R + Km * (Kb + Kp) );
plot(Kp, K); hold on;

% unit DC Gain for perfect steady-state tracking
plot(Kp, ones(size(Kp)), '--k'); 

xlabel('Gain of Proportional Controller, Kp');
ylabel('DC Gain, K');
title('P-Controller Gain influence over the DC Gain');

ylim([0 1.2])

grid on;