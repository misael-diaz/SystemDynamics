% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Course: System Dynamics and Control                       ME 3030-21 WI19
% Author: Prof. M Diaz-Maldonado                          Date: Dec/28/2019
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Synopsis:
% Example 10-6-A:
% Closed-loop response of mass-damper system controlled via a Proportional
% controller. Code shows that the response is "slow" if the controlled
% system is critically damped or overdamped (small Kp values). In the next
% example it will be shown that the response of its underdamped counterpart
% is much faster.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
clear;
close all;
clc;

% reference signal, step-input position, m 
A = 0.1;

% plot the reference step input
figure(1);
t = linspace(0, 200, 1e4);
% step input, u(t) = x_ref(t) = 0.1 U(t)
u = A * (t>=0);
plot(t, u, '--k'); hold on;

% P-Controller Gain (for Kp = 0.01125 the system is critically damped)
Kp = [0.005, 0.01, 0.01125];

for n=1:2 % (only the first two values correspond to an overdamped system)
    
    % second-order dynamic system constants: 
    %                y'' + a1 y' + a0 y = b0 u(t) 
    %
    a1 = 0.3;
    a0 = 2*Kp(n);
    % non-homogeneous constant
    b0 = 2*Kp(n);
    
    % natural frequency, rad/s
    omega_n = sqrt(a0);
    % damping ratio
    DR = a1 / (2*sqrt(a0));
    
    if (DR>1)
        % plots the response only if the system is overdamped
        % the conditional statement is not strictly needed, however
        % it has been added to stress that for other values of Kp the
        % system becomes underdamped and that means that y(t) is
        % described by an underdamped second-order response to a step input
        
        % roots of characteristic equation obtained by quadratic equation
        r1 = 0.5 * (-a1 + sqrt(a1^2 - 4 * a0));
        r2 = 0.5 * (-a1 - sqrt(a1^2 - 4 * a0));
        
        % closed-loop response (overdamped system) for a system that starts
        % at rest (y=y'=0 at t=0)
        y = (b0 * A / a0) * ( 1 + (r2/r1) / ( 1 - r2/r1 ) *...
            exp(r1*t) - 1/(1-r2/r1) * exp(r2*t) );
        plot(t, y); hold on; 
    end
 
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% critically damped case (Kp = 0.01125 V/m) 
r1 = -0.5 * a1;
% complete response, y(t), for a critically damped system which starts from
% rest (y=y'=0 at t=0).
y = (b0 * A/ a0) * (1 - exp(r1*t) + r1 * t .* exp(r1*t));
plot(t, y);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% final figure options
xlabel('Time, sec');
ylabel('Closed-loop Response, y(t) = x(t), m');
title('Response P-Controlled Mass-Damper System');
grid on;

ylim([0, 0.12]);
legend('step-input, x_{ref}(t)', 'K_p = 0.005 V/m (overdamped)', ... 
    'K_p = 0.01 V/m (overdamped)',...
    'K_p = 0.01125 V/m \bf{(critically damped)}','location', 'SE');