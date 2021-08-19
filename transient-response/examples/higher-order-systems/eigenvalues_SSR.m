% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Course: System Dynamics and Control                       ME 3030-21 WI19
% Author: Prof. M Diaz-Maldonado                          Date: Dec/27/2019
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Synopsis:
% Shows how to use the eig() built-in function to obtain the eigenvalues of
% the system matrix of the SSR of a linear dynamic system. It's shown that
% the eigenvalues are identical to the roots of the characteristic
% equation.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
clear;
close all;
clc;

% get the roots from the characteristic equation for the eigenvalues
C=[1 3 8 12];
r=roots(C)

% get the eigenvalues by invoking the built-in fun eig() on the system
% matrix
A = [ 0  1  0
      0  0  1
    -12 -8 -3];
r=eig(A) % eigenvalues are identical to the roots of the characteristic equation