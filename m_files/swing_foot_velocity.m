function [vh,vv] = swing_foot_velocity(x)
% SWING_FOOT_VELOCITY
%   [VH,VV] = SWING_FOOT_VELOCITY(X)  

%Eric Westervelt
%2016 Version: Peter Minh
%05-Dec-2016 17:03:23

modelP;

[nn,mm] = size(x);

if nn ~= 1 & mm ~= 1
  vh = x(:,6).*(L3.*cos(x(:,1) + x(:,5)) + L4.*cos(x(:,1) + x(:,3) + x(:,5))) - x(:,7).*(L3.*cos(x(:,2) + x(:,5)) + L4.*cos(x(:,2) + x(:,4) + x(:,5))) + x(:,10).*(L3.*cos(x(:,1) + x(:,5)) - L3.*cos(x(:,2) + x(:,5)) + L4.*cos(x(:,1) + x(:,3) + x(:,5)) - L4.*cos(x(:,2) + x(:,4) + x(:,5))) + L4.*x(:,8).*cos(x(:,1) + x(:,3) + x(:,5)) - L4.*x(:,9).*cos(x(:,2) + x(:,4) + x(:,5));
  vv = x(:,6).*(L3.*sin(x(:,1) + x(:,5)) + L4.*sin(x(:,1) + x(:,3) + x(:,5))) - x(:,7).*(L3.*sin(x(:,2) + x(:,5)) + L4.*sin(x(:,2) + x(:,4) + x(:,5))) + x(:,10).*(L3.*sin(x(:,1) + x(:,5)) - L3.*sin(x(:,2) + x(:,5)) + L4.*sin(x(:,1) + x(:,3) + x(:,5)) - L4.*sin(x(:,2) + x(:,4) + x(:,5))) + L4.*x(:,8).*sin(x(:,1) + x(:,3) + x(:,5)) - L4.*x(:,9).*sin(x(:,2) + x(:,4) + x(:,5));
else
  vh = x(6)*(L3*cos(x(1) + x(5)) + L4*cos(x(1) + x(3) + x(5))) - x(7)*(L3*cos(x(2) + x(5)) + L4*cos(x(2) + x(4) + x(5))) + x(10)*(L3*cos(x(1) + x(5)) - L3*cos(x(2) + x(5)) + L4*cos(x(1) + x(3) + x(5)) - L4*cos(x(2) + x(4) + x(5))) + L4*x(8)*cos(x(1) + x(3) + x(5)) - L4*x(9)*cos(x(2) + x(4) + x(5));
  vv = x(6)*(L3*sin(x(1) + x(5)) + L4*sin(x(1) + x(3) + x(5))) - x(7)*(L3*sin(x(2) + x(5)) + L4*sin(x(2) + x(4) + x(5))) + x(10)*(L3*sin(x(1) + x(5)) - L3*sin(x(2) + x(5)) + L4*sin(x(1) + x(3) + x(5)) - L4*sin(x(2) + x(4) + x(5))) + L4*x(8)*sin(x(1) + x(3) + x(5)) - L4*x(9)*sin(x(2) + x(4) + x(5));
end