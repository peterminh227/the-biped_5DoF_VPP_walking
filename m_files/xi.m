function [out] = xi(in)
% XI
%   [OUT] = XI(IN)  

%Eric Westervelt
%2016 Version: Peter Minh
%05-Dec-2016 17:03:24

xi = zeros(5);
xi(1,2) = 1;
xi(2,1) = 1;
xi(3,4) = 1;
xi(4,3) = 1;
xi(5,5) = 1;
out = xi*in;