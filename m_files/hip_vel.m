function [vHh,vHv] = hip_vel(x)
% HIP_VEL
%   [VHH,VHV] = HIP_VEL(X)  

%Eric Westervelt
%2016 Version: Peter Minh
%05-Dec-2016 17:03:23

modelP;

[nn,mm] = size(x);

if nn ~= 1 & mm ~= 1
  vHh = x(:,10).*(L3.*cos(x(:,1) + x(:,5)) + L4.*cos(x(:,1) + x(:,3) + x(:,5))) + x(:,6).*(L3.*cos(x(:,1) + x(:,5)) + L4.*cos(x(:,1) + x(:,3) + x(:,5))) + L4.*x(:,8).*cos(x(:,1) + x(:,3) + x(:,5));
  vHv = x(:,10).*(L3.*sin(x(:,1) + x(:,5)) + L4.*sin(x(:,1) + x(:,3) + x(:,5))) + x(:,6).*(L3.*sin(x(:,1) + x(:,5)) + L4.*sin(x(:,1) + x(:,3) + x(:,5))) + L4.*x(:,8).*sin(x(:,1) + x(:,3) + x(:,5));
else
  vHh = x(10)*(L3*cos(x(1) + x(5)) + L4*cos(x(1) + x(3) + x(5))) + x(6)*(L3*cos(x(1) + x(5)) + L4*cos(x(1) + x(3) + x(5))) + L4*x(8)*cos(x(1) + x(3) + x(5));
  vHv = x(10)*(L3*sin(x(1) + x(5)) + L4*sin(x(1) + x(3) + x(5))) + x(6)*(L3*sin(x(1) + x(5)) + L4*sin(x(1) + x(3) + x(5))) + L4*x(8)*sin(x(1) + x(3) + x(5));
end
