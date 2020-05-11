function [alpha_times_beta] = alpha_beta(z1,a,b,m)
% ALPHA_BETA
%   [ALPHA_TIMES_BETA] = ALPHA_BETA(Z1,A,B,M)  

%Eric Westervelt
%11-Nov-2016 14:55:10

alpha_val = alpha(z1,a,b,m);
beta_val = beta(z1,a,b,m);

alpha_times_beta = alpha_val.*beta_val;