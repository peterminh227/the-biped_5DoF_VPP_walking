function [h] = head_height(x)
% HEAD_HEIGHT
%   [H] = HEAD_HEIGHT(X)  

%Eric Westervelt
%2016 Version: Peter Minh
%05-Dec-2016 17:03:23

modelP;

[nn,mm] = size(x);

if nn ~= 1 & mm ~= 1
  h = L1.*cos(x(:,5)) - L3.*cos(x(:,1) + x(:,5)) - L4.*cos(x(:,1) + x(:,3) + x(:,5));
else
  h = L1*cos(x(5)) - L3*cos(x(1) + x(5)) - L4*cos(x(1) + x(3) + x(5));
end
