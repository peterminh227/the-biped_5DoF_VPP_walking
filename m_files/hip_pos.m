function [pHh,pHv] = hip_pos(x)
% HIP_POS
%   [PHH,PHV] = HIP_POS(X)  

%Eric Westervelt
%2016 Version: Peter Minh
%05-Dec-2016 17:03:23

modelP;

[nn,mm] = size(x);

if nn ~= 1 & mm ~= 1
  pHh = L3.*sin(x(:,1) + x(:,5)) + L4.*sin(x(:,1) + x(:,3) + x(:,5));
  pHv = - L3.*cos(x(:,1) + x(:,5)) - L4.*cos(x(:,1) + x(:,3) + x(:,5));
else
  pHh = L3*sin(x(1) + x(5)) + L4*sin(x(1) + x(3) + x(5));
  pHv = - L3*cos(x(1) + x(5)) - L4*cos(x(1) + x(3) + x(5));
end
