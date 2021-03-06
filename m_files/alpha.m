function [alpha] = alpha(z1,a,b,m)
% ALPHA
%   [ALPHA] = ALPHA(Z1,A,B,M)  

%Eric Westervelt
%11-Nov-2016 14:55:10

modelP;

a0=a(1); a1=a(2); a2=a(3); a3=a(4); a4=a(5); a5=a(6); a6=a(7); 
b0=a(8); b1=a(9); b2=a(10); b3=a(11); b4=a(12); b5=a(13); b6=a(14); 
c0=a(15); c1=a(16); c2=a(17); c3=a(18); c4=a(19); c5=a(20); c6=a(21); 
d0=a(22); d1=a(23); d2=a(24); d3=a(25); d4=a(26); d5=a(27); d6=a(28); 

hd1 = a0.*((b - z1)./m + 1).^6 + (a6.*(b - z1).^6)./m.^6 + (15.*a2.* ...
      ((b - z1)./m + 1).^4.*(b - z1).^2)./m.^2 - (20.*a3.*((b - ...
       z1)./m + 1).^3.*(b - z1).^3)./m.^3 + (15.*a4.*((b - z1)./m + ...
       1).^2.*(b - z1).^4)./m.^4 - (6.*a1.*((b - z1)./m + 1).^5.*(b - ...
       z1))./m - (6.*a5.*((b - z1)./m + 1).*(b - z1).^5)./m.^5;
hd2 = b0.*((b - z1)./m + 1).^6 + (b6.*(b - z1).^6)./m.^6 + (15.*b2.* ...
      ((b - z1)./m + 1).^4.*(b - z1).^2)./m.^2 - (20.*b3.*((b - ...
       z1)./m + 1).^3.*(b - z1).^3)./m.^3 + (15.*b4.*((b - z1)./m + ...
       1).^2.*(b - z1).^4)./m.^4 - (6.*b1.*((b - z1)./m + 1).^5.*(b - ...
       z1))./m - (6.*b5.*((b - z1)./m + 1).*(b - z1).^5)./m.^5;
hd3 = c0.*((b - z1)./m + 1).^6 + (c6.*(b - z1).^6)./m.^6 + (15.*c2.* ...
      ((b - z1)./m + 1).^4.*(b - z1).^2)./m.^2 - (20.*c3.*((b - ...
       z1)./m + 1).^3.*(b - z1).^3)./m.^3 + (15.*c4.*((b - z1)./m + ...
       1).^2.*(b - z1).^4)./m.^4 - (6.*c1.*((b - z1)./m + 1).^5.*(b - ...
       z1))./m - (6.*c5.*((b - z1)./m + 1).*(b - z1).^5)./m.^5;
hd4 = d0.*((b - z1)./m + 1).^6 + (d6.*(b - z1).^6)./m.^6 + (15.*d2.* ...
      ((b - z1)./m + 1).^4.*(b - z1).^2)./m.^2 - (20.*d3.*((b - ...
       z1)./m + 1).^3.*(b - z1).^3)./m.^3 + (15.*d4.*((b - z1)./m + ...
       1).^2.*(b - z1).^4)./m.^4 - (6.*d1.*((b - z1)./m + 1).^5.*(b - ...
       z1))./m - (6.*d5.*((b - z1)./m + 1).*(b - z1).^5)./m.^5;

jac_hd1 = (6.*a1.*((b - z1)./m + 1).^5)./m - (6.*a0.*((b - z1)./m + ...
           1).^5)./m + (6.*a5.*(b - z1).^5)./m.^6 - (6.*a6.*(b - ...
           z1).^5)./m.^6 - (15.*a2.*((b - z1)./m + 1).^4.*(2.*b - 2.* ...
          z1))./m.^2 - (60.*a2.*((b - z1)./m + 1).^3.*(b - z1).^2)./m.^3 + ...
           (60.*a3.*((b - z1)./m + 1).^3.*(b - z1).^2)./m.^3 + (60.*a3.* ...
          ((b - z1)./m + 1).^2.*(b - z1).^3)./m.^4 - (60.*a4.*((b - ...
           z1)./m + 1).^2.*(b - z1).^3)./m.^4 + (30.*a1.*((b - z1)./m + ...
           1).^4.*(b - z1))./m.^2 - (30.*a4.*((b - z1)./m + 1).*(b - ...
           z1).^4)./m.^5 + (30.*a5.*((b - z1)./m + 1).*(b - z1).^4)./m.^5;
jac_hd2 = (6.*b1.*((b - z1)./m + 1).^5)./m - (6.*b0.*((b - z1)./m + ...
           1).^5)./m + (6.*b5.*(b - z1).^5)./m.^6 - (6.*b6.*(b - ...
           z1).^5)./m.^6 - (15.*b2.*((b - z1)./m + 1).^4.*(2.*b - 2.* ...
          z1))./m.^2 - (60.*b2.*((b - z1)./m + 1).^3.*(b - z1).^2)./m.^3 + ...
           (60.*b3.*((b - z1)./m + 1).^3.*(b - z1).^2)./m.^3 + (60.*b3.* ...
          ((b - z1)./m + 1).^2.*(b - z1).^3)./m.^4 - (60.*b4.*((b - ...
           z1)./m + 1).^2.*(b - z1).^3)./m.^4 + (30.*b1.*((b - z1)./m + ...
           1).^4.*(b - z1))./m.^2 - (30.*b4.*((b - z1)./m + 1).*(b - ...
           z1).^4)./m.^5 + (30.*b5.*((b - z1)./m + 1).*(b - z1).^4)./m.^5;
jac_hd3 = (6.*c1.*((b - z1)./m + 1).^5)./m - (6.*c0.*((b - z1)./m + ...
           1).^5)./m + (6.*c5.*(b - z1).^5)./m.^6 - (6.*c6.*(b - ...
           z1).^5)./m.^6 - (15.*c2.*((b - z1)./m + 1).^4.*(2.*b - 2.* ...
          z1))./m.^2 - (60.*c2.*((b - z1)./m + 1).^3.*(b - z1).^2)./m.^3 + ...
           (60.*c3.*((b - z1)./m + 1).^3.*(b - z1).^2)./m.^3 + (60.*c3.* ...
          ((b - z1)./m + 1).^2.*(b - z1).^3)./m.^4 - (60.*c4.*((b - ...
           z1)./m + 1).^2.*(b - z1).^3)./m.^4 + (30.*c1.*((b - z1)./m + ...
           1).^4.*(b - z1))./m.^2 - (30.*c4.*((b - z1)./m + 1).*(b - ...
           z1).^4)./m.^5 + (30.*c5.*((b - z1)./m + 1).*(b - z1).^4)./m.^5;
jac_hd4 = (6.*d1.*((b - z1)./m + 1).^5)./m - (6.*d0.*((b - z1)./m + ...
           1).^5)./m + (6.*d5.*(b - z1).^5)./m.^6 - (6.*d6.*(b - ...
           z1).^5)./m.^6 - (15.*d2.*((b - z1)./m + 1).^4.*(2.*b - 2.* ...
          z1))./m.^2 - (60.*d2.*((b - z1)./m + 1).^3.*(b - z1).^2)./m.^3 + ...
           (60.*d3.*((b - z1)./m + 1).^3.*(b - z1).^2)./m.^3 + (60.*d3.* ...
          ((b - z1)./m + 1).^2.*(b - z1).^3)./m.^4 - (60.*d4.*((b - ...
           z1)./m + 1).^2.*(b - z1).^3)./m.^4 + (30.*d1.*((b - z1)./m + ...
           1).^4.*(b - z1))./m.^2 - (30.*d4.*((b - z1)./m + 1).*(b - ...
           z1).^4)./m.^5 + (30.*d5.*((b - z1)./m + 1).*(b - z1).^4)./m.^5;

alpha = XX3.*jac_hd2 - 2.*XX3 - 2.*XX4 - XX1.*jac_hd1 - (XX1.* ...
        jac_hd3)./2 - XX3.*jac_hd1 - XX1 - XX4.*jac_hd1 - XX3.*jac_hd3 + ...
         XX4.*jac_hd2 + XX4.*jac_hd4 - L3.^2.*M1 - L4.^2.*M1 - 2.*L3.^2.* ...
        M3 - 2.*L3.^2.*M4 - 2.*L4.^2.*M3 - 2.*L4.^2.*M4 - M1.*MY1.^2 - ...
         2.*M3.*MY3.^2 - 2.*M4.*MY4.^2 - M1.*MZ1.^2 - 2.*M3.*MZ3.^2 - 2.* ...
        M4.*MZ4.^2 + 2.*L3.*M3.*MZ3 + 2.*L4.*M4.*MZ4 + 2.*L3.^2.*M4.* ...
        cos(hd1 - hd2) - (L3.^2.*M1.*jac_hd3)./2 - L3.^2.*M4.*jac_hd1 + ...
         (L4.^2.*M1.*jac_hd3)./2 - L3.^2.*M3.*jac_hd3 + L3.^2.*M4.* ...
        jac_hd2 - L3.^2.*M4.*jac_hd3 + L4.^2.*M3.*jac_hd3 + L4.^2.*M4.* ...
        jac_hd3 - M1.*MY1.^2.*jac_hd1 - (M1.*MY1.^2.*jac_hd3)./2 - M3.* ...
        MY3.^2.*jac_hd1 + M3.*MY3.^2.*jac_hd2 - M3.*MY3.^2.*jac_hd3 - ...
         M4.*MY4.^2.*jac_hd1 + M4.*MY4.^2.*jac_hd2 + M4.*MY4.^2.* ...
        jac_hd4 - M1.*MZ1.^2.*jac_hd1 - (M1.*MZ1.^2.*jac_hd3)./2 - M3.* ...
        MZ3.^2.*jac_hd1 + M3.*MZ3.^2.*jac_hd2 - M3.*MZ3.^2.*jac_hd3 - ...
         M4.*MZ4.^2.*jac_hd1 + M4.*MZ4.^2.*jac_hd2 + M4.*MZ4.^2.* ...
        jac_hd4 + 2.*L3.*M4.*MZ4.*cos(hd1 - hd2 - hd4) - 2.*L3.*L4.*M1.* ...
        cos(hd3) - 4.*L3.*L4.*M3.*cos(hd3) - 2.*L3.*L4.*M4.*cos(hd3) - ...
         2.*L3.*M4.*MY4.*sin(hd1 - hd2 - hd4) + 2.*L3.*M1.*MZ1.* ...
        cos(hd1) + 2.*L4.*M3.*MZ3.*cos(hd3) - 2.*L3.*M4.*MZ4.*cos(hd4) - ...
         2.*L3.*M1.*MY1.*sin(hd1) - 2.*L4.*M3.*MY3.*sin(hd3) - 2.*L3.* ...
        M4.*MY4.*sin(hd4) + 2.*L3.*M3.*MZ3.*cos(hd1 - hd2) - 2.*L3.*M3.* ...
        MY3.*sin(hd1 - hd2) + 2.*L4.*M4.*MZ4.*cos(hd1 - hd2 + hd3 - ...
         hd4) - 2.*L4.*M4.*MY4.*sin(hd1 - hd2 + hd3 - hd4) + 2.*L3.*L4.* ...
        M4.*cos(hd1 - hd2 + hd3) + 2.*L4.*M3.*MZ3.*cos(hd1 - hd2 + hd3) + ...
         L3.*M3.*MZ3.*jac_hd3 - L4.*M4.*MZ4.*jac_hd3 - 2.*L4.*M3.*MY3.* ...
        sin(hd1 - hd2 + hd3) + L3.^2.*M4.*jac_hd1.*cos(hd1 - hd2) - ...
         L3.^2.*M4.*jac_hd2.*cos(hd1 - hd2) + L3.^2.*M4.*jac_hd3.* ...
        cos(hd1 - hd2) + 2.*L4.*M1.*MZ1.*cos(hd1 + hd3) - 2.*L4.*M1.* ...
        MY1.*sin(hd1 + hd3) + L3.*L4.*M4.*jac_hd1.*cos(hd1 - hd2 + hd3) - ...
         L3.*L4.*M4.*jac_hd2.*cos(hd1 - hd2 + hd3) + L4.*M3.*MZ3.* ...
        jac_hd1.*cos(hd1 - hd2 + hd3) - L4.*M3.*MZ3.*jac_hd2.*cos(hd1 - ...
         hd2 + hd3) - L4.*M3.*MY3.*jac_hd1.*sin(hd1 - hd2 + hd3) + L4.* ...
        M3.*MY3.*jac_hd2.*sin(hd1 - hd2 + hd3) + L4.*M1.*MZ1.*jac_hd1.* ...
        cos(hd1 + hd3) - L4.*M1.*MY1.*jac_hd1.*sin(hd1 + hd3) + L3.*M4.* ...
        MZ4.*jac_hd1.*cos(hd1 - hd2 - hd4) - L3.*M4.*MZ4.*jac_hd2.* ...
        cos(hd1 - hd2 - hd4) + L3.*M4.*MZ4.*jac_hd3.*cos(hd1 - hd2 - ...
         hd4) - L3.*M4.*MZ4.*jac_hd4.*cos(hd1 - hd2 - hd4) - L3.*M4.* ...
        MY4.*jac_hd1.*sin(hd1 - hd2 - hd4) + L3.*M4.*MY4.*jac_hd2.* ...
        sin(hd1 - hd2 - hd4) - L3.*M4.*MY4.*jac_hd3.*sin(hd1 - hd2 - ...
         hd4) + L3.*M4.*MY4.*jac_hd4.*sin(hd1 - hd2 - hd4) + L3.*M1.* ...
        MZ1.*jac_hd1.*cos(hd1) + L3.*M1.*MZ1.*jac_hd3.*cos(hd1) - 2.*L3.* ...
        M4.*MZ4.*jac_hd1.*cos(hd4) + 2.*L3.*M4.*MZ4.*jac_hd2.*cos(hd4) - ...
         L3.*M4.*MZ4.*jac_hd3.*cos(hd4) + L3.*M4.*MZ4.*jac_hd4.* ...
        cos(hd4) - L3.*M1.*MY1.*jac_hd1.*sin(hd1) - L3.*M1.*MY1.* ...
        jac_hd3.*sin(hd1) - 2.*L3.*M4.*MY4.*jac_hd1.*sin(hd4) + 2.*L3.* ...
        M4.*MY4.*jac_hd2.*sin(hd4) - L3.*M4.*MY4.*jac_hd3.*sin(hd4) + ...
         L3.*M4.*MY4.*jac_hd4.*sin(hd4) + L3.*M3.*MZ3.*jac_hd1.*cos(hd1 - ...
         hd2) - L3.*M3.*MZ3.*jac_hd2.*cos(hd1 - hd2) + L3.*M3.*MZ3.* ...
        jac_hd3.*cos(hd1 - hd2) - L3.*M3.*MY3.*jac_hd1.*sin(hd1 - hd2) + ...
         L3.*M3.*MY3.*jac_hd2.*sin(hd1 - hd2) - L3.*M3.*MY3.*jac_hd3.* ...
        sin(hd1 - hd2) + L4.*M4.*MZ4.*jac_hd1.*cos(hd1 - hd2 + hd3 - ...
         hd4) - L4.*M4.*MZ4.*jac_hd2.*cos(hd1 - hd2 + hd3 - hd4) - L4.* ...
        M4.*MZ4.*jac_hd4.*cos(hd1 - hd2 + hd3 - hd4) - L4.*M4.*MY4.* ...
        jac_hd1.*sin(hd1 - hd2 + hd3 - hd4) + L4.*M4.*MY4.*jac_hd2.* ...
        sin(hd1 - hd2 + hd3 - hd4) + L4.*M4.*MY4.*jac_hd4.*sin(hd1 - ...
         hd2 + hd3 - hd4);