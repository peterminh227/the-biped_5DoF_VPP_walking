function [pCh,pCv] = CoM_pos(x)
% COM_POS
%   [PCH,PCV] = COM_POS(X)  

%Eric Westervelt
%2016 Version: Peter Minh
%05-Dec-2016 17:03:23

modelP;

[nn,mm] = size(x);

if nn ~= 1 & mm ~= 1
  pCh = (M4.*(L3.*sin(x(:,1) + x(:,5)) - L3.*sin(x(:,2) + x(:,5)) + MY4.*cos(x(:,2) + x(:,4) + x(:,5)) + L4.*sin(x(:,1) + x(:,3) + x(:,5)) - MZ4.*sin(x(:,2) + x(:,4) + x(:,5))) + M1.*(L3.*sin(x(:,1) + x(:,5)) + MY1.*cos(x(:,5)) - MZ1.*sin(x(:,5)) + L4.*sin(x(:,1) + x(:,3) + x(:,5))) + M4.*(MY4.*cos(x(:,1) + x(:,3) + x(:,5)) + L4.*sin(x(:,1) + x(:,3) + x(:,5)) - MZ4.*sin(x(:,1) + x(:,3) + x(:,5))) + M3.*(MY3.*cos(x(:,1) + x(:,5)) + L3.*sin(x(:,1) + x(:,5)) - MZ3.*sin(x(:,1) + x(:,5)) + L4.*sin(x(:,1) + x(:,3) + x(:,5))) + M3.*(MY3.*cos(x(:,2) + x(:,5)) + L3.*sin(x(:,1) + x(:,5)) - MZ3.*sin(x(:,2) + x(:,5)) + L4.*sin(x(:,1) + x(:,3) + x(:,5))))./(M1 + 2.*M3 + 2.*M4);
  pCv = -(M1.*(L3.*cos(x(:,1) + x(:,5)) - MZ1.*cos(x(:,5)) - MY1.*sin(x(:,5)) + L4.*cos(x(:,1) + x(:,3) + x(:,5))) - M4.*(L3.*cos(x(:,2) + x(:,5)) - L3.*cos(x(:,1) + x(:,5)) - L4.*cos(x(:,1) + x(:,3) + x(:,5)) + MZ4.*cos(x(:,2) + x(:,4) + x(:,5)) + MY4.*sin(x(:,2) + x(:,4) + x(:,5))) - M4.*(MZ4.*cos(x(:,1) + x(:,3) + x(:,5)) - L4.*cos(x(:,1) + x(:,3) + x(:,5)) + MY4.*sin(x(:,1) + x(:,3) + x(:,5))) + M3.*(L3.*cos(x(:,1) + x(:,5)) - MZ3.*cos(x(:,1) + x(:,5)) - MY3.*sin(x(:,1) + x(:,5)) + L4.*cos(x(:,1) + x(:,3) + x(:,5))) + M3.*(L3.*cos(x(:,1) + x(:,5)) - MZ3.*cos(x(:,2) + x(:,5)) - MY3.*sin(x(:,2) + x(:,5)) + L4.*cos(x(:,1) + x(:,3) + x(:,5))))./(M1 + 2.*M3 + 2.*M4);
else
  pCh = (M4*(L3*sin(x(1) + x(5)) - L3*sin(x(2) + x(5)) + MY4*cos(x(2) + x(4) + x(5)) + L4*sin(x(1) + x(3) + x(5)) - MZ4*sin(x(2) + x(4) + x(5))) + M1*(L3*sin(x(1) + x(5)) + MY1*cos(x(5)) - MZ1*sin(x(5)) + L4*sin(x(1) + x(3) + x(5))) + M4*(MY4*cos(x(1) + x(3) + x(5)) + L4*sin(x(1) + x(3) + x(5)) - MZ4*sin(x(1) + x(3) + x(5))) + M3*(MY3*cos(x(1) + x(5)) + L3*sin(x(1) + x(5)) - MZ3*sin(x(1) + x(5)) + L4*sin(x(1) + x(3) + x(5))) + M3*(MY3*cos(x(2) + x(5)) + L3*sin(x(1) + x(5)) - MZ3*sin(x(2) + x(5)) + L4*sin(x(1) + x(3) + x(5))))/(M1 + 2*M3 + 2*M4);
  pCv = -(M1*(L3*cos(x(1) + x(5)) - MZ1*cos(x(5)) - MY1*sin(x(5)) + L4*cos(x(1) + x(3) + x(5))) - M4*(L3*cos(x(2) + x(5)) - L3*cos(x(1) + x(5)) - L4*cos(x(1) + x(3) + x(5)) + MZ4*cos(x(2) + x(4) + x(5)) + MY4*sin(x(2) + x(4) + x(5))) - M4*(MZ4*cos(x(1) + x(3) + x(5)) - L4*cos(x(1) + x(3) + x(5)) + MY4*sin(x(1) + x(3) + x(5))) + M3*(L3*cos(x(1) + x(5)) - MZ3*cos(x(1) + x(5)) - MY3*sin(x(1) + x(5)) + L4*cos(x(1) + x(3) + x(5))) + M3*(L3*cos(x(1) + x(5)) - MZ3*cos(x(2) + x(5)) - MY3*sin(x(2) + x(5)) + L4*cos(x(1) + x(3) + x(5))))/(M1 + 2*M3 + 2*M4);
end
