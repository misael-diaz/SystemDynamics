% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Course: System Dynamics and Control                       ME 3030-21 WI19
% Author: Prof. M Diaz-Maldonado                          Date: Jan/15/2020
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Synopsis:
% Constructs the Bode Plot for a thrid-order mass-spring-damper system.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
clear;
close all;
clc;

% 
k = 1; b = 5; m = 1;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% transfer functions
G0 = @(s) 1./(m*s.^2+k);
G1 = @(s) k./(b*s+k);
G2 = @(s) k./(m*s.^2+k);

G = @(s) G0(s) ./ (1-G1(s).*G2(s));

% ampplitude ratio
K = @(s) abs(G(s));
% phase angle
phi = @(s) angle(G(s));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Produce the Bode Plot
figure(1);

omega = logspace(-1,2,500);
subplot(2,1,1);
semilogx(omega, 20*log10(K(j*omega)) );
ylabel('Magnitude of G(jw), dB');
xlabel('Frequency, \omega, rad/s');
title('Third-order mass-spring-damper system');
%ylim([-50, 10])
grid on;

subplot(2,1,2);
semilogx(omega, 180 * phi(j*omega) / pi);
ylabel('Phase of G(jw), deg');
xlabel('Frequency, \omega, rad/s');
%ylim([-95, 5]);

grid on;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
 