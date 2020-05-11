function [h] = step_length(x)
% STEP_LENGTH
%   [H] = STEP_LENGTH(X)  

%Eric Westervelt
%2016 Version: Peter Minh
%05-Dec-2016 17:03:23

modelP;

q31L=x(1); q32L=x(2); q41L=x(3); q42L=x(4); q1L=x(5);
h = L3*sin(q31L + q1L) - L3*sin(q32L + q1L) + L4*sin(q31L + q41L + ...
 q1L) - L4*sin(q32L + q42L + q1L);