% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Course: System Dynamics and Control                       ME 3030-21 WI19
% Author: Prof. M Diaz-Maldonado                          Date: Dec/28/2019
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Synopsis:
% Example 10-6-A-2:
% Closed-loop response of a underdamped P-controlled mass-damper system.
% Although the response reaches steady state faster (~30 seconds) than its
% overdamped counterpart, the response exhibits an undesired overshoot as
% the controller gain, Kp, is increased. It's noted that the settling time
% remains constant for increasing values Kp.
% 
% Code reproduces Figure 10.19 of Kluever's textbook.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
clear;
close all;
clc;

% reference signal, step-input position, m 
A = 0.1;

% plot the reference step input
figure(1);
t = linspace(0, 30, 1e4);
% step input, u(t) = x_ref(t) = 0.1 U(t)
u = A * (t>=0);
plot(t, u, '--k'); hold on;

% P-Controller Gain 
Kp = [0.2, 1, 5];

for n=1:length(Kp)
    
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
    
    if (DR>0 || DR < 1)
        % plots the response, y(t) for a underdamped system
        
        
        % real part of the characteristic roots
        alpha = -(DR * omega_n);
        % imaginary part of the characteristic roots
        beta = omega_n * sqrt(1-DR^2);
        
        % closed-loop response, y(t)
        y = (b0 * A / a0) * ( 1 - exp(alpha * t) .* ( cos(beta*t) -...
            alpha/beta * sin(beta*t) ) );
        plot(t, y); hold on;
    end
 
end

% final figure options
xlabel('Time, sec');
ylabel('Closed-loop Response, y(t) = x(t), m');
title('Response Underdamped P-Controlled Mass-Damper System');
grid on;

ylim([0, 0.2]);
legend('step-input, x_{ref}(t)', 'K_p = 0.2 V/m', 'K_p = 1.0 V/m',...
    'K_p = 5.0 V/m');