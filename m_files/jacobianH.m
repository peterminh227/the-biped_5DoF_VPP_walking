function [jac_h] = jacobianH(q,dq,a,b,m)
% JACOBIANH
%   [JAC_H] = JACOBIANH(Q,DQ,A,B,M)  

%Eric Westervelt
%12-Nov-2016 22:00:19

q31L=q(1);q32L=q(2);q41L=q(3);q42L=q(4);q1L=q(5);
dq31L=dq(1);dq32L=dq(2);dq41L=dq(3);dq42L=dq(4);dq1L=dq(5);

a0=a(1); a1=a(2); a2=a(3); a3=a(4); a4=a(5); a5=a(6); a6=a(7); 
b0=a(8); b1=a(9); b2=a(10); b3=a(11); b4=a(12); b5=a(13); b6=a(14); 
c0=a(15); c1=a(16); c2=a(17); c3=a(18); c4=a(19); c5=a(20); c6=a(21); 
d0=a(22); d1=a(23); d2=a(24); d3=a(25); d4=a(26); d5=a(27); d6=a(28); 


s = (- q31L - q41L/2 - q1L - b)/m;
if s > 1
  s = 1;
elseif s < 0
  s = 0;
end
z1 = m*s + b;

dz1 = - dq1L - dq31L - dq41L/2;

jac_h = zeros(4,5);

jac_h(1,1) = (6*a1*((b - z1)/m + 1)^5)/m - (6*a0*((b - z1)/m + 1)^5)/m + ...
              (6*a5*(b - z1)^5)/m^6 - (6*a6*(b - z1)^5)/m^6 - (15*a2* ...
             ((b - z1)/m + 1)^4*(2*b - 2*z1))/m^2 - (60*a2*((b - z1)/m + ...
              1)^3*(b - z1)^2)/m^3 + (60*a3*((b - z1)/m + 1)^3*(b - ...
              z1)^2)/m^3 + (60*a3*((b - z1)/m + 1)^2*(b - z1)^3)/m^4 - ...
              (60*a4*((b - z1)/m + 1)^2*(b - z1)^3)/m^4 + (30*a1*((b - ...
              z1)/m + 1)^4*(b - z1))/m^2 - (30*a4*((b - z1)/m + 1)*(b - ...
              z1)^4)/m^5 + (30*a5*((b - z1)/m + 1)*(b - z1)^4)/m^5 + 1;
jac_h(1,3) = (3*a1*((b - z1)/m + 1)^5)/m - (3*a0*((b - z1)/m + 1)^5)/m + ...
              (3*a5*(b - z1)^5)/m^6 - (3*a6*(b - z1)^5)/m^6 - (15*a2* ...
             ((b - z1)/m + 1)^4*(2*b - 2*z1))/(2*m^2) - (30*a2*((b - ...
              z1)/m + 1)^3*(b - z1)^2)/m^3 + (30*a3*((b - z1)/m + 1)^3* ...
             (b - z1)^2)/m^3 + (30*a3*((b - z1)/m + 1)^2*(b - ...
              z1)^3)/m^4 - (30*a4*((b - z1)/m + 1)^2*(b - z1)^3)/m^4 + ...
              (15*a1*((b - z1)/m + 1)^4*(b - z1))/m^2 - (15*a4*((b - ...
              z1)/m + 1)*(b - z1)^4)/m^5 + (15*a5*((b - z1)/m + 1)*(b - ...
              z1)^4)/m^5;
jac_h(1,5) = (6*a1*((b - z1)/m + 1)^5)/m - (6*a0*((b - z1)/m + 1)^5)/m + ...
              (6*a5*(b - z1)^5)/m^6 - (6*a6*(b - z1)^5)/m^6 - (15*a2* ...
             ((b - z1)/m + 1)^4*(2*b - 2*z1))/m^2 - (60*a2*((b - z1)/m + ...
              1)^3*(b - z1)^2)/m^3 + (60*a3*((b - z1)/m + 1)^3*(b - ...
              z1)^2)/m^3 + (60*a3*((b - z1)/m + 1)^2*(b - z1)^3)/m^4 - ...
              (60*a4*((b - z1)/m + 1)^2*(b - z1)^3)/m^4 + (30*a1*((b - ...
              z1)/m + 1)^4*(b - z1))/m^2 - (30*a4*((b - z1)/m + 1)*(b - ...
              z1)^4)/m^5 + (30*a5*((b - z1)/m + 1)*(b - z1)^4)/m^5;
jac_h(2,1) = (6*b1*((b - z1)/m + 1)^5)/m - (6*b0*((b - z1)/m + 1)^5)/m + ...
              (6*b5*(b - z1)^5)/m^6 - (6*b6*(b - z1)^5)/m^6 - (15*b2* ...
             ((b - z1)/m + 1)^4*(2*b - 2*z1))/m^2 - (60*b2*((b - z1)/m + ...
              1)^3*(b - z1)^2)/m^3 + (60*b3*((b - z1)/m + 1)^3*(b - ...
              z1)^2)/m^3 + (60*b3*((b - z1)/m + 1)^2*(b - z1)^3)/m^4 - ...
              (60*b4*((b - z1)/m + 1)^2*(b - z1)^3)/m^4 + (30*b1*((b - ...
              z1)/m + 1)^4*(b - z1))/m^2 - (30*b4*((b - z1)/m + 1)*(b - ...
              z1)^4)/m^5 + (30*b5*((b - z1)/m + 1)*(b - z1)^4)/m^5;
jac_h(2,2) = 1;
jac_h(2,3) = (3*b1*((b - z1)/m + 1)^5)/m - (3*b0*((b - z1)/m + 1)^5)/m + ...
              (3*b5*(b - z1)^5)/m^6 - (3*b6*(b - z1)^5)/m^6 - (15*b2* ...
             ((b - z1)/m + 1)^4*(2*b - 2*z1))/(2*m^2) - (30*b2*((b - ...
              z1)/m + 1)^3*(b - z1)^2)/m^3 + (30*b3*((b - z1)/m + 1)^3* ...
             (b - z1)^2)/m^3 + (30*b3*((b - z1)/m + 1)^2*(b - ...
              z1)^3)/m^4 - (30*b4*((b - z1)/m + 1)^2*(b - z1)^3)/m^4 + ...
              (15*b1*((b - z1)/m + 1)^4*(b - z1))/m^2 - (15*b4*((b - ...
              z1)/m + 1)*(b - z1)^4)/m^5 + (15*b5*((b - z1)/m + 1)*(b - ...
              z1)^4)/m^5;
jac_h(2,5) = (6*b1*((b - z1)/m + 1)^5)/m - (6*b0*((b - z1)/m + 1)^5)/m + ...
              (6*b5*(b - z1)^5)/m^6 - (6*b6*(b - z1)^5)/m^6 - (15*b2* ...
             ((b - z1)/m + 1)^4*(2*b - 2*z1))/m^2 - (60*b2*((b - z1)/m + ...
              1)^3*(b - z1)^2)/m^3 + (60*b3*((b - z1)/m + 1)^3*(b - ...
              z1)^2)/m^3 + (60*b3*((b - z1)/m + 1)^2*(b - z1)^3)/m^4 - ...
              (60*b4*((b - z1)/m + 1)^2*(b - z1)^3)/m^4 + (30*b1*((b - ...
              z1)/m + 1)^4*(b - z1))/m^2 - (30*b4*((b - z1)/m + 1)*(b - ...
              z1)^4)/m^5 + (30*b5*((b - z1)/m + 1)*(b - z1)^4)/m^5;
jac_h(3,1) = (6*c1*((b - z1)/m + 1)^5)/m - (6*c0*((b - z1)/m + 1)^5)/m + ...
              (6*c5*(b - z1)^5)/m^6 - (6*c6*(b - z1)^5)/m^6 - (15*c2* ...
             ((b - z1)/m + 1)^4*(2*b - 2*z1))/m^2 - (60*c2*((b - z1)/m + ...
              1)^3*(b - z1)^2)/m^3 + (60*c3*((b - z1)/m + 1)^3*(b - ...
              z1)^2)/m^3 + (60*c3*((b - z1)/m + 1)^2*(b - z1)^3)/m^4 - ...
              (60*c4*((b - z1)/m + 1)^2*(b - z1)^3)/m^4 + (30*c1*((b - ...
              z1)/m + 1)^4*(b - z1))/m^2 - (30*c4*((b - z1)/m + 1)*(b - ...
              z1)^4)/m^5 + (30*c5*((b - z1)/m + 1)*(b - z1)^4)/m^5;
jac_h(3,3) = (3*c1*((b - z1)/m + 1)^5)/m - (3*c0*((b - z1)/m + 1)^5)/m + ...
              (3*c5*(b - z1)^5)/m^6 - (3*c6*(b - z1)^5)/m^6 - (15*c2* ...
             ((b - z1)/m + 1)^4*(2*b - 2*z1))/(2*m^2) - (30*c2*((b - ...
              z1)/m + 1)^3*(b - z1)^2)/m^3 + (30*c3*((b - z1)/m + 1)^3* ...
             (b - z1)^2)/m^3 + (30*c3*((b - z1)/m + 1)^2*(b - ...
              z1)^3)/m^4 - (30*c4*((b - z1)/m + 1)^2*(b - z1)^3)/m^4 + ...
              (15*c1*((b - z1)/m + 1)^4*(b - z1))/m^2 - (15*c4*((b - ...
              z1)/m + 1)*(b - z1)^4)/m^5 + (15*c5*((b - z1)/m + 1)*(b - ...
              z1)^4)/m^5 + 1;
jac_h(3,5) = (6*c1*((b - z1)/m + 1)^5)/m - (6*c0*((b - z1)/m + 1)^5)/m + ...
              (6*c5*(b - z1)^5)/m^6 - (6*c6*(b - z1)^5)/m^6 - (15*c2* ...
             ((b - z1)/m + 1)^4*(2*b - 2*z1))/m^2 - (60*c2*((b - z1)/m + ...
              1)^3*(b - z1)^2)/m^3 + (60*c3*((b - z1)/m + 1)^3*(b - ...
              z1)^2)/m^3 + (60*c3*((b - z1)/m + 1)^2*(b - z1)^3)/m^4 - ...
              (60*c4*((b - z1)/m + 1)^2*(b - z1)^3)/m^4 + (30*c1*((b - ...
              z1)/m + 1)^4*(b - z1))/m^2 - (30*c4*((b - z1)/m + 1)*(b - ...
              z1)^4)/m^5 + (30*c5*((b - z1)/m + 1)*(b - z1)^4)/m^5;
jac_h(4,1) = (6*d1*((b - z1)/m + 1)^5)/m - (6*d0*((b - z1)/m + 1)^5)/m + ...
              (6*d5*(b - z1)^5)/m^6 - (6*d6*(b - z1)^5)/m^6 - (15*d2* ...
             ((b - z1)/m + 1)^4*(2*b - 2*z1))/m^2 - (60*d2*((b - z1)/m + ...
              1)^3*(b - z1)^2)/m^3 + (60*d3*((b - z1)/m + 1)^3*(b - ...
              z1)^2)/m^3 + (60*d3*((b - z1)/m + 1)^2*(b - z1)^3)/m^4 - ...
              (60*d4*((b - z1)/m + 1)^2*(b - z1)^3)/m^4 + (30*d1*((b - ...
              z1)/m + 1)^4*(b - z1))/m^2 - (30*d4*((b - z1)/m + 1)*(b - ...
              z1)^4)/m^5 + (30*d5*((b - z1)/m + 1)*(b - z1)^4)/m^5;
jac_h(4,3) = (3*d1*((b - z1)/m + 1)^5)/m - (3*d0*((b - z1)/m + 1)^5)/m + ...
              (3*d5*(b - z1)^5)/m^6 - (3*d6*(b - z1)^5)/m^6 - (15*d2* ...
             ((b - z1)/m + 1)^4*(2*b - 2*z1))/(2*m^2) - (30*d2*((b - ...
              z1)/m + 1)^3*(b - z1)^2)/m^3 + (30*d3*((b - z1)/m + 1)^3* ...
             (b - z1)^2)/m^3 + (30*d3*((b - z1)/m + 1)^2*(b - ...
              z1)^3)/m^4 - (30*d4*((b - z1)/m + 1)^2*(b - z1)^3)/m^4 + ...
              (15*d1*((b - z1)/m + 1)^4*(b - z1))/m^2 - (15*d4*((b - ...
              z1)/m + 1)*(b - z1)^4)/m^5 + (15*d5*((b - z1)/m + 1)*(b - ...
              z1)^4)/m^5;
jac_h(4,4) = 1;
jac_h(4,5) = (6*d1*((b - z1)/m + 1)^5)/m - (6*d0*((b - z1)/m + 1)^5)/m + ...
              (6*d5*(b - z1)^5)/m^6 - (6*d6*(b - z1)^5)/m^6 - (15*d2* ...
             ((b - z1)/m + 1)^4*(2*b - 2*z1))/m^2 - (60*d2*((b - z1)/m + ...
              1)^3*(b - z1)^2)/m^3 + (60*d3*((b - z1)/m + 1)^3*(b - ...
              z1)^2)/m^3 + (60*d3*((b - z1)/m + 1)^2*(b - z1)^3)/m^4 - ...
              (60*d4*((b - z1)/m + 1)^2*(b - z1)^3)/m^4 + (30*d1*((b - ...
              z1)/m + 1)^4*(b - z1))/m^2 - (30*d4*((b - z1)/m + 1)*(b - ...
              z1)^4)/m^5 + (30*d5*((b - z1)/m + 1)*(b - z1)^4)/m^5;