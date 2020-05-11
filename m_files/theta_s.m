function [z1] = theta_s(x)
% THETA_S
%   [Z1] = THETA_S(X)  

%Eric Westervelt
%12-Nov-2016 22:00:20

[nn,mm] = size(x);

if nn ~= 1 & mm ~= 1
  z1 = - x(:,1) - x(:,3)./2 - x(:,5);
else
  z1 = - x(1) - x(3)/2 - x(5);
end
