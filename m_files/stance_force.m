function [f_tan,f_norm] = stance_force(x,ddq)
% STANCE_FORCE
%   [F_TAN,F_NORM] = STANCE_FORCE(X,DDQ)  

%Eric Westervelt
%2016 Version: Peter Minh
%16-Nov-2016 10:40:29

modelP;

% tangential force
f_tan = L3*M1*ddq(5)*cos(x(1) + x(5)) + 2*L3*M3*ddq(5)*cos(x(1) + x(5)) + ...
         L3*M4*ddq(5)*cos(x(1) + x(5)) - L3*M4*ddq(5)*cos(x(2) + x(5)) + ...
         L3*M1*ddq(1)*cos(x(1) + x(5)) + 2*L3*M3*ddq(1)*cos(x(1) + ...
         x(5)) + L3*M4*ddq(1)*cos(x(1) + x(5)) - L3*M4*ddq(2)*cos(x(2) + ...
         x(5)) - M3*MZ3*ddq(5)*cos(x(1) + x(5)) - M3*MZ3*ddq(5)* ...
        cos(x(2) + x(5)) - M3*MZ3*ddq(1)*cos(x(1) + x(5)) - M3*MZ3* ...
        ddq(2)*cos(x(2) + x(5)) - M3*MY3*ddq(5)*sin(x(1) + x(5)) - M3* ...
        MY3*ddq(5)*sin(x(2) + x(5)) - M3*MY3*ddq(1)*sin(x(1) + x(5)) - ...
         M3*MY3*ddq(2)*sin(x(2) + x(5)) - M1*MZ1*ddq(5)*cos(x(5)) - M1* ...
        MY1*ddq(5)*sin(x(5)) - M3*MY3*x(10)^2*cos(x(1) + x(5)) - M3*MY3* ...
        x(10)^2*cos(x(2) + x(5)) - M3*MY3*x(6)^2*cos(x(1) + x(5)) - M3* ...
        MY3*x(7)^2*cos(x(2) + x(5)) - L3*M1*x(10)^2*sin(x(1) + x(5)) - 2* ...
        L3*M3*x(10)^2*sin(x(1) + x(5)) - L3*M4*x(10)^2*sin(x(1) + x(5)) + ...
         L3*M4*x(10)^2*sin(x(2) + x(5)) - L3*M1*x(6)^2*sin(x(1) + x(5)) - ...
         2*L3*M3*x(6)^2*sin(x(1) + x(5)) - L3*M4*x(6)^2*sin(x(1) + ...
         x(5)) + L3*M4*x(7)^2*sin(x(2) + x(5)) + M3*MZ3*x(10)^2* ...
        sin(x(1) + x(5)) + M3*MZ3*x(10)^2*sin(x(2) + x(5)) + M3*MZ3* ...
        x(6)^2*sin(x(1) + x(5)) + M3*MZ3*x(7)^2*sin(x(2) + x(5)) + L4*M1* ...
        ddq(5)*cos(x(1) + x(3) + x(5)) + 2*L4*M3*ddq(5)*cos(x(1) + x(3) + ...
         x(5)) + 2*L4*M4*ddq(5)*cos(x(1) + x(3) + x(5)) + L4*M1*ddq(1)* ...
        cos(x(1) + x(3) + x(5)) + 2*L4*M3*ddq(1)*cos(x(1) + x(3) + ...
         x(5)) + 2*L4*M4*ddq(1)*cos(x(1) + x(3) + x(5)) + L4*M1*ddq(3)* ...
        cos(x(1) + x(3) + x(5)) + 2*L4*M3*ddq(3)*cos(x(1) + x(3) + ...
         x(5)) + 2*L4*M4*ddq(3)*cos(x(1) + x(3) + x(5)) - M4*MZ4*ddq(5)* ...
        cos(x(1) + x(3) + x(5)) - M4*MZ4*ddq(5)*cos(x(2) + x(4) + x(5)) - ...
         M4*MZ4*ddq(1)*cos(x(1) + x(3) + x(5)) - M4*MZ4*ddq(2)*cos(x(2) + ...
         x(4) + x(5)) - M4*MZ4*ddq(3)*cos(x(1) + x(3) + x(5)) - M4*MZ4* ...
        ddq(4)*cos(x(2) + x(4) + x(5)) - M4*MY4*ddq(5)*sin(x(1) + x(3) + ...
         x(5)) - M4*MY4*ddq(5)*sin(x(2) + x(4) + x(5)) - M4*MY4*ddq(1)* ...
        sin(x(1) + x(3) + x(5)) - M4*MY4*ddq(2)*sin(x(2) + x(4) + x(5)) - ...
         M4*MY4*ddq(3)*sin(x(1) + x(3) + x(5)) - M4*MY4*ddq(4)*sin(x(2) + ...
         x(4) + x(5)) - M1*MY1*x(10)^2*cos(x(5)) + M1*MZ1*x(10)^2* ...
        sin(x(5)) - M4*MY4*x(10)^2*cos(x(1) + x(3) + x(5)) - M4*MY4* ...
        x(10)^2*cos(x(2) + x(4) + x(5)) - M4*MY4*x(6)^2*cos(x(1) + x(3) + ...
         x(5)) - M4*MY4*x(7)^2*cos(x(2) + x(4) + x(5)) - M4*MY4*x(8)^2* ...
        cos(x(1) + x(3) + x(5)) - M4*MY4*x(9)^2*cos(x(2) + x(4) + x(5)) - ...
         L4*M1*x(10)^2*sin(x(1) + x(3) + x(5)) - 2*L4*M3*x(10)^2* ...
        sin(x(1) + x(3) + x(5)) - 2*L4*M4*x(10)^2*sin(x(1) + x(3) + ...
         x(5)) - L4*M1*x(6)^2*sin(x(1) + x(3) + x(5)) - 2*L4*M3*x(6)^2* ...
        sin(x(1) + x(3) + x(5)) - 2*L4*M4*x(6)^2*sin(x(1) + x(3) + ...
         x(5)) - L4*M1*x(8)^2*sin(x(1) + x(3) + x(5)) - 2*L4*M3*x(8)^2* ...
        sin(x(1) + x(3) + x(5)) - 2*L4*M4*x(8)^2*sin(x(1) + x(3) + ...
         x(5)) + M4*MZ4*x(10)^2*sin(x(1) + x(3) + x(5)) + M4*MZ4*x(10)^2* ...
        sin(x(2) + x(4) + x(5)) + M4*MZ4*x(6)^2*sin(x(1) + x(3) + x(5)) + ...
         M4*MZ4*x(7)^2*sin(x(2) + x(4) + x(5)) + M4*MZ4*x(8)^2*sin(x(1) + ...
         x(3) + x(5)) + M4*MZ4*x(9)^2*sin(x(2) + x(4) + x(5)) - 2*M3*MY3* ...
        x(10)*x(6)*cos(x(1) + x(5)) - 2*M3*MY3*x(10)*x(7)*cos(x(2) + ...
         x(5)) - 2*L3*M1*x(10)*x(6)*sin(x(1) + x(5)) - 4*L3*M3*x(10)* ...
        x(6)*sin(x(1) + x(5)) - 2*L3*M4*x(10)*x(6)*sin(x(1) + x(5)) + 2* ...
        L3*M4*x(10)*x(7)*sin(x(2) + x(5)) + 2*M3*MZ3*x(10)*x(6)* ...
        sin(x(1) + x(5)) + 2*M3*MZ3*x(10)*x(7)*sin(x(2) + x(5)) - 2*M4* ...
        MY4*x(10)*x(6)*cos(x(1) + x(3) + x(5)) - 2*M4*MY4*x(10)*x(7)* ...
        cos(x(2) + x(4) + x(5)) - 2*M4*MY4*x(10)*x(8)*cos(x(1) + x(3) + ...
         x(5)) - 2*M4*MY4*x(10)*x(9)*cos(x(2) + x(4) + x(5)) - 2*M4*MY4* ...
        x(6)*x(8)*cos(x(1) + x(3) + x(5)) - 2*M4*MY4*x(7)*x(9)*cos(x(2) + ...
         x(4) + x(5)) - 2*L4*M1*x(10)*x(6)*sin(x(1) + x(3) + x(5)) - 4* ...
        L4*M3*x(10)*x(6)*sin(x(1) + x(3) + x(5)) - 4*L4*M4*x(10)*x(6)* ...
        sin(x(1) + x(3) + x(5)) - 2*L4*M1*x(10)*x(8)*sin(x(1) + x(3) + ...
         x(5)) - 4*L4*M3*x(10)*x(8)*sin(x(1) + x(3) + x(5)) - 4*L4*M4* ...
        x(10)*x(8)*sin(x(1) + x(3) + x(5)) - 2*L4*M1*x(6)*x(8)*sin(x(1) + ...
         x(3) + x(5)) - 4*L4*M3*x(6)*x(8)*sin(x(1) + x(3) + x(5)) - 4*L4* ...
        M4*x(6)*x(8)*sin(x(1) + x(3) + x(5)) + 2*M4*MZ4*x(10)*x(6)* ...
        sin(x(1) + x(3) + x(5)) + 2*M4*MZ4*x(10)*x(7)*sin(x(2) + x(4) + ...
         x(5)) + 2*M4*MZ4*x(10)*x(8)*sin(x(1) + x(3) + x(5)) + 2*M4*MZ4* ...
        x(10)*x(9)*sin(x(2) + x(4) + x(5)) + 2*M4*MZ4*x(6)*x(8)* ...
        sin(x(1) + x(3) + x(5)) + 2*M4*MZ4*x(7)*x(9)*sin(x(2) + x(4) + ...
         x(5));

% normal force
f_norm = M1*g + 2*M3*g + 2*M4*g + M3*MY3*ddq(5)*cos(x(1) + x(5)) + M3*MY3* ...
         ddq(5)*cos(x(2) + x(5)) + M3*MY3*ddq(1)*cos(x(1) + x(5)) + M3* ...
         MY3*ddq(2)*cos(x(2) + x(5)) + L3*M1*ddq(5)*sin(x(1) + x(5)) + 2* ...
         L3*M3*ddq(5)*sin(x(1) + x(5)) + L3*M4*ddq(5)*sin(x(1) + x(5)) - ...
          L3*M4*ddq(5)*sin(x(2) + x(5)) + L3*M1*ddq(1)*sin(x(1) + x(5)) + ...
          2*L3*M3*ddq(1)*sin(x(1) + x(5)) + L3*M4*ddq(1)*sin(x(1) + ...
          x(5)) - L3*M4*ddq(2)*sin(x(2) + x(5)) - M3*MZ3*ddq(5)*sin(x(1) + ...
          x(5)) - M3*MZ3*ddq(5)*sin(x(2) + x(5)) - M3*MZ3*ddq(1)* ...
         sin(x(1) + x(5)) - M3*MZ3*ddq(2)*sin(x(2) + x(5)) + M1*MY1* ...
         ddq(5)*cos(x(5)) - M1*MZ1*ddq(5)*sin(x(5)) + L3*M1*x(10)^2* ...
         cos(x(1) + x(5)) + 2*L3*M3*x(10)^2*cos(x(1) + x(5)) + L3*M4* ...
         x(10)^2*cos(x(1) + x(5)) - L3*M4*x(10)^2*cos(x(2) + x(5)) + L3* ...
         M1*x(6)^2*cos(x(1) + x(5)) + 2*L3*M3*x(6)^2*cos(x(1) + x(5)) + ...
          L3*M4*x(6)^2*cos(x(1) + x(5)) - L3*M4*x(7)^2*cos(x(2) + x(5)) - ...
          M3*MZ3*x(10)^2*cos(x(1) + x(5)) - M3*MZ3*x(10)^2*cos(x(2) + ...
          x(5)) - M3*MZ3*x(6)^2*cos(x(1) + x(5)) - M3*MZ3*x(7)^2* ...
         cos(x(2) + x(5)) - M3*MY3*x(10)^2*sin(x(1) + x(5)) - M3*MY3* ...
         x(10)^2*sin(x(2) + x(5)) - M3*MY3*x(6)^2*sin(x(1) + x(5)) - M3* ...
         MY3*x(7)^2*sin(x(2) + x(5)) + M4*MY4*ddq(5)*cos(x(1) + x(3) + ...
          x(5)) + M4*MY4*ddq(5)*cos(x(2) + x(4) + x(5)) + M4*MY4*ddq(1)* ...
         cos(x(1) + x(3) + x(5)) + M4*MY4*ddq(2)*cos(x(2) + x(4) + x(5)) + ...
          M4*MY4*ddq(3)*cos(x(1) + x(3) + x(5)) + M4*MY4*ddq(4)*cos(x(2) + ...
          x(4) + x(5)) + L4*M1*ddq(5)*sin(x(1) + x(3) + x(5)) + 2*L4*M3* ...
         ddq(5)*sin(x(1) + x(3) + x(5)) + 2*L4*M4*ddq(5)*sin(x(1) + x(3) + ...
          x(5)) + L4*M1*ddq(1)*sin(x(1) + x(3) + x(5)) + 2*L4*M3*ddq(1)* ...
         sin(x(1) + x(3) + x(5)) + 2*L4*M4*ddq(1)*sin(x(1) + x(3) + ...
          x(5)) + L4*M1*ddq(3)*sin(x(1) + x(3) + x(5)) + 2*L4*M3*ddq(3)* ...
         sin(x(1) + x(3) + x(5)) + 2*L4*M4*ddq(3)*sin(x(1) + x(3) + ...
          x(5)) - M4*MZ4*ddq(5)*sin(x(1) + x(3) + x(5)) - M4*MZ4*ddq(5)* ...
         sin(x(2) + x(4) + x(5)) - M4*MZ4*ddq(1)*sin(x(1) + x(3) + x(5)) - ...
          M4*MZ4*ddq(2)*sin(x(2) + x(4) + x(5)) - M4*MZ4*ddq(3)*sin(x(1) + ...
          x(3) + x(5)) - M4*MZ4*ddq(4)*sin(x(2) + x(4) + x(5)) - M1*MZ1* ...
         x(10)^2*cos(x(5)) - M1*MY1*x(10)^2*sin(x(5)) + L4*M1*x(10)^2* ...
         cos(x(1) + x(3) + x(5)) + 2*L4*M3*x(10)^2*cos(x(1) + x(3) + ...
          x(5)) + 2*L4*M4*x(10)^2*cos(x(1) + x(3) + x(5)) + L4*M1*x(6)^2* ...
         cos(x(1) + x(3) + x(5)) + 2*L4*M3*x(6)^2*cos(x(1) + x(3) + ...
          x(5)) + 2*L4*M4*x(6)^2*cos(x(1) + x(3) + x(5)) + L4*M1*x(8)^2* ...
         cos(x(1) + x(3) + x(5)) + 2*L4*M3*x(8)^2*cos(x(1) + x(3) + ...
          x(5)) + 2*L4*M4*x(8)^2*cos(x(1) + x(3) + x(5)) - M4*MZ4*x(10)^2* ...
         cos(x(1) + x(3) + x(5)) - M4*MZ4*x(10)^2*cos(x(2) + x(4) + ...
          x(5)) - M4*MZ4*x(6)^2*cos(x(1) + x(3) + x(5)) - M4*MZ4*x(7)^2* ...
         cos(x(2) + x(4) + x(5)) - M4*MZ4*x(8)^2*cos(x(1) + x(3) + x(5)) - ...
          M4*MZ4*x(9)^2*cos(x(2) + x(4) + x(5)) - M4*MY4*x(10)^2* ...
         sin(x(1) + x(3) + x(5)) - M4*MY4*x(10)^2*sin(x(2) + x(4) + ...
          x(5)) - M4*MY4*x(6)^2*sin(x(1) + x(3) + x(5)) - M4*MY4*x(7)^2* ...
         sin(x(2) + x(4) + x(5)) - M4*MY4*x(8)^2*sin(x(1) + x(3) + x(5)) - ...
          M4*MY4*x(9)^2*sin(x(2) + x(4) + x(5)) + 2*L3*M1*x(10)*x(6)* ...
         cos(x(1) + x(5)) + 4*L3*M3*x(10)*x(6)*cos(x(1) + x(5)) + 2*L3*M4* ...
         x(10)*x(6)*cos(x(1) + x(5)) - 2*L3*M4*x(10)*x(7)*cos(x(2) + ...
          x(5)) - 2*M3*MZ3*x(10)*x(6)*cos(x(1) + x(5)) - 2*M3*MZ3*x(10)* ...
         x(7)*cos(x(2) + x(5)) - 2*M3*MY3*x(10)*x(6)*sin(x(1) + x(5)) - 2* ...
         M3*MY3*x(10)*x(7)*sin(x(2) + x(5)) + 2*L4*M1*x(10)*x(6)* ...
         cos(x(1) + x(3) + x(5)) + 4*L4*M3*x(10)*x(6)*cos(x(1) + x(3) + ...
          x(5)) + 4*L4*M4*x(10)*x(6)*cos(x(1) + x(3) + x(5)) + 2*L4*M1* ...
         x(10)*x(8)*cos(x(1) + x(3) + x(5)) + 4*L4*M3*x(10)*x(8)* ...
         cos(x(1) + x(3) + x(5)) + 4*L4*M4*x(10)*x(8)*cos(x(1) + x(3) + ...
          x(5)) + 2*L4*M1*x(6)*x(8)*cos(x(1) + x(3) + x(5)) + 4*L4*M3* ...
         x(6)*x(8)*cos(x(1) + x(3) + x(5)) + 4*L4*M4*x(6)*x(8)*cos(x(1) + ...
          x(3) + x(5)) - 2*M4*MZ4*x(10)*x(6)*cos(x(1) + x(3) + x(5)) - 2* ...
         M4*MZ4*x(10)*x(7)*cos(x(2) + x(4) + x(5)) - 2*M4*MZ4*x(10)*x(8)* ...
         cos(x(1) + x(3) + x(5)) - 2*M4*MZ4*x(10)*x(9)*cos(x(2) + x(4) + ...
          x(5)) - 2*M4*MZ4*x(6)*x(8)*cos(x(1) + x(3) + x(5)) - 2*M4*MZ4* ...
         x(7)*x(9)*cos(x(2) + x(4) + x(5)) - 2*M4*MY4*x(10)*x(6)* ...
         sin(x(1) + x(3) + x(5)) - 2*M4*MY4*x(10)*x(7)*sin(x(2) + x(4) + ...
          x(5)) - 2*M4*MY4*x(10)*x(8)*sin(x(1) + x(3) + x(5)) - 2*M4*MY4* ...
         x(10)*x(9)*sin(x(2) + x(4) + x(5)) - 2*M4*MY4*x(6)*x(8)* ...
         sin(x(1) + x(3) + x(5)) - 2*M4*MY4*x(7)*x(9)*sin(x(2) + x(4) + ...
          x(5));