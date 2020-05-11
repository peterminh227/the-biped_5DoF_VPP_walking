function [dz1] = dtheta_s(x)
% DTHETA_S
%   [DZ1] = DTHETA_S(X)  

%Eric Westervelt
%12-Nov-2016 22:00:20

[nn,mm] = size(x);

if nn ~= 1 & mm ~= 1
  dz1 = - x(:,10) - x(:,6) - x(:,8)./2;
else
  dz1 = - x(10) - x(6) - x(8)/2;
end
