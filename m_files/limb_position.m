function [out] = limb_position(x,pH_horiz)
% LIMB_POSITION
%   [OUT] = LIMB_POSITION(X,PH_HORIZ)  

%Eric Westervelt
%2016 Version: Peter Minh
%05-Dec-2016 17:03:23

modelP;

[nn,mm] = size(x);
if nn == 10, k = mm; else, k = nn; end

if nn ~= 1 & mm ~= 1
  out.pH1 = pH_horiz*ones(k,1);
  out.pH2 = - L3.*cos(x(:,1) + x(:,5)) - L4.*cos(x(:,1) + x(:,3) + x(:,5));
  out.pT1 = out.pH1 + MY1.*cos(x(:,5)) - MZ1.*sin(x(:,5));
  out.pT2 = out.pH2 + MZ1.*cos(x(:,5)) + MY1.*sin(x(:,5));
  out.pCoM1 = out.pH1 + (M4.*(L3.*sin(x(:,1) + x(:,5)) - L3.*sin(x(:,2) + x(:,5)) + MY4.*cos(x(:,2) + x(:,4) + x(:,5)) + L4.*sin(x(:,1) + x(:,3) + x(:,5)) - MZ4.*sin(x(:,2) + x(:,4) + x(:,5))) + M1.*(L3.*sin(x(:,1) + x(:,5)) + MY1.*cos(x(:,5)) - MZ1.*sin(x(:,5)) + L4.*sin(x(:,1) + x(:,3) + x(:,5))) + M4.*(MY4.*cos(x(:,1) + x(:,3) + x(:,5)) + L4.*sin(x(:,1) + x(:,3) + x(:,5)) - MZ4.*sin(x(:,1) + x(:,3) + x(:,5))) + M3.*(MY3.*cos(x(:,1) + x(:,5)) + L3.*sin(x(:,1) + x(:,5)) - MZ3.*sin(x(:,1) + x(:,5)) + L4.*sin(x(:,1) + x(:,3) + x(:,5))) + M3.*(MY3.*cos(x(:,2) + x(:,5)) + L3.*sin(x(:,1) + x(:,5)) - MZ3.*sin(x(:,2) + x(:,5)) + L4.*sin(x(:,1) + x(:,3) + x(:,5))))./(M1 + 2.*M3 + 2.*M4) - L4.*sin(x(:,1) + x(:,3) + x(:,5)) - L3.*sin(x(:,1) + x(:,5));
  out.pCoM2 = out.pH2 + L3.*cos(x(:,1) + x(:,5)) - (M1.*(L3.*cos(x(:,1) + x(:,5)) - MZ1.*cos(x(:,5)) - MY1.*sin(x(:,5)) + L4.*cos(x(:,1) + x(:,3) + x(:,5))) - M4.*(L3.*cos(x(:,2) + x(:,5)) - L3.*cos(x(:,1) + x(:,5)) - L4.*cos(x(:,1) + x(:,3) + x(:,5)) + MZ4.*cos(x(:,2) + x(:,4) + x(:,5)) + MY4.*sin(x(:,2) + x(:,4) + x(:,5))) - M4.*(MZ4.*cos(x(:,1) + x(:,3) + x(:,5)) - L4.*cos(x(:,1) + x(:,3) + x(:,5)) + MY4.*sin(x(:,1) + x(:,3) + x(:,5))) + M3.*(L3.*cos(x(:,1) + x(:,5)) - MZ3.*cos(x(:,1) + x(:,5)) - MY3.*sin(x(:,1) + x(:,5)) + L4.*cos(x(:,1) + x(:,3) + x(:,5))) + M3.*(L3.*cos(x(:,1) + x(:,5)) - MZ3.*cos(x(:,2) + x(:,5)) - MY3.*sin(x(:,2) + x(:,5)) + L4.*cos(x(:,1) + x(:,3) + x(:,5))))./(M1 + 2.*M3 + 2.*M4) + L4.*cos(x(:,1) + x(:,3) + x(:,5));
  out.pVPP1 = out.pH1 + (M4.*(L3.*sin(x(:,1) + x(:,5)) - L3.*sin(x(:,2) + x(:,5)) + MY4.*cos(x(:,2) + x(:,4) + x(:,5)) + L4.*sin(x(:,1) + x(:,3) + x(:,5)) - MZ4.*sin(x(:,2) + x(:,4) + x(:,5))) + M1.*(L3.*sin(x(:,1) + x(:,5)) + MY1.*cos(x(:,5)) - MZ1.*sin(x(:,5)) + L4.*sin(x(:,1) + x(:,3) + x(:,5))) + M4.*(MY4.*cos(x(:,1) + x(:,3) + x(:,5)) + L4.*sin(x(:,1) + x(:,3) + x(:,5)) - MZ4.*sin(x(:,1) + x(:,3) + x(:,5))) + M3.*(MY3.*cos(x(:,1) + x(:,5)) + L3.*sin(x(:,1) + x(:,5)) - MZ3.*sin(x(:,1) + x(:,5)) + L4.*sin(x(:,1) + x(:,3) + x(:,5))) + M3.*(MY3.*cos(x(:,2) + x(:,5)) + L3.*sin(x(:,1) + x(:,5)) - MZ3.*sin(x(:,2) + x(:,5)) + L4.*sin(x(:,1) + x(:,3) + x(:,5))))./(M1 + 2.*M3 + 2.*M4) - rVPP.*sin(gam) - L4.*sin(x(:,1) + x(:,3) + x(:,5)) - L3.*sin(x(:,1) + x(:,5));
  out.pVPP2 = out.pH2 + L3.*cos(x(:,1) + x(:,5)) - (M1.*(L3.*cos(x(:,1) + x(:,5)) - MZ1.*cos(x(:,5)) - MY1.*sin(x(:,5)) + L4.*cos(x(:,1) + x(:,3) + x(:,5))) - M4.*(L3.*cos(x(:,2) + x(:,5)) - L3.*cos(x(:,1) + x(:,5)) - L4.*cos(x(:,1) + x(:,3) + x(:,5)) + MZ4.*cos(x(:,2) + x(:,4) + x(:,5)) + MY4.*sin(x(:,2) + x(:,4) + x(:,5))) - M4.*(MZ4.*cos(x(:,1) + x(:,3) + x(:,5)) - L4.*cos(x(:,1) + x(:,3) + x(:,5)) + MY4.*sin(x(:,1) + x(:,3) + x(:,5))) + M3.*(L3.*cos(x(:,1) + x(:,5)) - MZ3.*cos(x(:,1) + x(:,5)) - MY3.*sin(x(:,1) + x(:,5)) + L4.*cos(x(:,1) + x(:,3) + x(:,5))) + M3.*(L3.*cos(x(:,1) + x(:,5)) - MZ3.*cos(x(:,2) + x(:,5)) - MY3.*sin(x(:,2) + x(:,5)) + L4.*cos(x(:,1) + x(:,3) + x(:,5))))./(M1 + 2.*M3 + 2.*M4) + rVPP.*cos(gam) + L4.*cos(x(:,1) + x(:,3) + x(:,5));
  out.pHead1 = out.pH1 + -L1.*sin(x(:,5));
  out.pHead2 = out.pH2 + L1.*cos(x(:,5));
  out.pFem11 = out.pH1 + MY3.*cos(x(:,1) + x(:,5)) - MZ3.*sin(x(:,1) + x(:,5));
  out.pFem12 = out.pH2 + MZ3.*cos(x(:,1) + x(:,5)) + MY3.*sin(x(:,1) + x(:,5));
  out.pFem21 = out.pH1 + MY3.*cos(x(:,2) + x(:,5)) - MZ3.*sin(x(:,2) + x(:,5));
  out.pFem22 = out.pH2 + MZ3.*cos(x(:,2) + x(:,5)) + MY3.*sin(x(:,2) + x(:,5));
  out.pG11 = out.pH1 + -L3.*sin(x(:,1) + x(:,5));
  out.pG12 = out.pH2 + L3.*cos(x(:,1) + x(:,5));
  out.pG21 = out.pH1 + -L3.*sin(x(:,2) + x(:,5));
  out.pG22 = out.pH2 + L3.*cos(x(:,2) + x(:,5));
  out.pTib11 = out.pH1 + MY4.*cos(x(:,1) + x(:,3) + x(:,5)) - L3.*sin(x(:,1) + x(:,5)) - MZ4.*sin(x(:,1) + x(:,3) + x(:,5));
  out.pTib12 = out.pH2 + L3.*cos(x(:,1) + x(:,5)) + MZ4.*cos(x(:,1) + x(:,3) + x(:,5)) + MY4.*sin(x(:,1) + x(:,3) + x(:,5));
  out.pTib21 = out.pH1 + MY4.*cos(x(:,2) + x(:,4) + x(:,5)) - L3.*sin(x(:,2) + x(:,5)) - MZ4.*sin(x(:,2) + x(:,4) + x(:,5));
  out.pTib22 = out.pH2 + L3.*cos(x(:,2) + x(:,5)) + MZ4.*cos(x(:,2) + x(:,4) + x(:,5)) + MY4.*sin(x(:,2) + x(:,4) + x(:,5));
  out.pFoot11 = out.pH1 + - L3.*sin(x(:,1) + x(:,5)) - L4.*sin(x(:,1) + x(:,3) + x(:,5));
  out.pFoot12 = out.pH2 + L3.*cos(x(:,1) + x(:,5)) + L4.*cos(x(:,1) + x(:,3) + x(:,5));
  out.pFoot21 = out.pH1 + - L3.*sin(x(:,2) + x(:,5)) - L4.*sin(x(:,2) + x(:,4) + x(:,5));
  out.pFoot22 = out.pH2 + L3.*cos(x(:,2) + x(:,5)) + L4.*cos(x(:,2) + x(:,4) + x(:,5));
else
  out.pH1 = pH_horiz;
  out.pH2 = - L3*cos(x(1) + x(5)) - L4*cos(x(1) + x(3) + x(5));
  out.pT1 = out.pH1 + MY1*cos(x(5)) - MZ1*sin(x(5));
  out.pT2 = out.pH2 + MZ1*cos(x(5)) + MY1*sin(x(5));
  out.pCoM1 = out.pH1 + (M4.*(L3.*sin(x(:,1) + x(:,5)) - L3.*sin(x(:,2) + x(:,5)) + MY4.*cos(x(:,2) + x(:,4) + x(:,5)) + L4.*sin(x(:,1) + x(:,3) + x(:,5)) - MZ4.*sin(x(:,2) + x(:,4) + x(:,5))) + M1.*(L3.*sin(x(:,1) + x(:,5)) + MY1.*cos(x(:,5)) - MZ1.*sin(x(:,5)) + L4.*sin(x(:,1) + x(:,3) + x(:,5))) + M4.*(MY4.*cos(x(:,1) + x(:,3) + x(:,5)) + L4.*sin(x(:,1) + x(:,3) + x(:,5)) - MZ4.*sin(x(:,1) + x(:,3) + x(:,5))) + M3.*(MY3.*cos(x(:,1) + x(:,5)) + L3.*sin(x(:,1) + x(:,5)) - MZ3.*sin(x(:,1) + x(:,5)) + L4.*sin(x(:,1) + x(:,3) + x(:,5))) + M3.*(MY3.*cos(x(:,2) + x(:,5)) + L3.*sin(x(:,1) + x(:,5)) - MZ3.*sin(x(:,2) + x(:,5)) + L4.*sin(x(:,1) + x(:,3) + x(:,5))))./(M1 + 2.*M3 + 2.*M4) - L4.*sin(x(:,1) + x(:,3) + x(:,5)) - L3.*sin(x(:,1) + x(:,5));
  out.pCoM2 = out.pH2 + L3.*cos(x(:,1) + x(:,5)) - (M1.*(L3.*cos(x(:,1) + x(:,5)) - MZ1.*cos(x(:,5)) - MY1.*sin(x(:,5)) + L4.*cos(x(:,1) + x(:,3) + x(:,5))) - M4.*(L3.*cos(x(:,2) + x(:,5)) - L3.*cos(x(:,1) + x(:,5)) - L4.*cos(x(:,1) + x(:,3) + x(:,5)) + MZ4.*cos(x(:,2) + x(:,4) + x(:,5)) + MY4.*sin(x(:,2) + x(:,4) + x(:,5))) - M4.*(MZ4.*cos(x(:,1) + x(:,3) + x(:,5)) - L4.*cos(x(:,1) + x(:,3) + x(:,5)) + MY4.*sin(x(:,1) + x(:,3) + x(:,5))) + M3.*(L3.*cos(x(:,1) + x(:,5)) - MZ3.*cos(x(:,1) + x(:,5)) - MY3.*sin(x(:,1) + x(:,5)) + L4.*cos(x(:,1) + x(:,3) + x(:,5))) + M3.*(L3.*cos(x(:,1) + x(:,5)) - MZ3.*cos(x(:,2) + x(:,5)) - MY3.*sin(x(:,2) + x(:,5)) + L4.*cos(x(:,1) + x(:,3) + x(:,5))))./(M1 + 2.*M3 + 2.*M4) + L4.*cos(x(:,1) + x(:,3) + x(:,5));
  out.pVPP1 = out.pH1 + (M4.*(L3.*sin(x(:,1) + x(:,5)) - L3.*sin(x(:,2) + x(:,5)) + MY4.*cos(x(:,2) + x(:,4) + x(:,5)) + L4.*sin(x(:,1) + x(:,3) + x(:,5)) - MZ4.*sin(x(:,2) + x(:,4) + x(:,5))) + M1.*(L3.*sin(x(:,1) + x(:,5)) + MY1.*cos(x(:,5)) - MZ1.*sin(x(:,5)) + L4.*sin(x(:,1) + x(:,3) + x(:,5))) + M4.*(MY4.*cos(x(:,1) + x(:,3) + x(:,5)) + L4.*sin(x(:,1) + x(:,3) + x(:,5)) - MZ4.*sin(x(:,1) + x(:,3) + x(:,5))) + M3.*(MY3.*cos(x(:,1) + x(:,5)) + L3.*sin(x(:,1) + x(:,5)) - MZ3.*sin(x(:,1) + x(:,5)) + L4.*sin(x(:,1) + x(:,3) + x(:,5))) + M3.*(MY3.*cos(x(:,2) + x(:,5)) + L3.*sin(x(:,1) + x(:,5)) - MZ3.*sin(x(:,2) + x(:,5)) + L4.*sin(x(:,1) + x(:,3) + x(:,5))))./(M1 + 2.*M3 + 2.*M4) - rVPP.*sin(gam) - L4.*sin(x(:,1) + x(:,3) + x(:,5)) - L3.*sin(x(:,1) + x(:,5));
  out.pVPP2 = out.pH2 + L3.*cos(x(:,1) + x(:,5)) - (M1.*(L3.*cos(x(:,1) + x(:,5)) - MZ1.*cos(x(:,5)) - MY1.*sin(x(:,5)) + L4.*cos(x(:,1) + x(:,3) + x(:,5))) - M4.*(L3.*cos(x(:,2) + x(:,5)) - L3.*cos(x(:,1) + x(:,5)) - L4.*cos(x(:,1) + x(:,3) + x(:,5)) + MZ4.*cos(x(:,2) + x(:,4) + x(:,5)) + MY4.*sin(x(:,2) + x(:,4) + x(:,5))) - M4.*(MZ4.*cos(x(:,1) + x(:,3) + x(:,5)) - L4.*cos(x(:,1) + x(:,3) + x(:,5)) + MY4.*sin(x(:,1) + x(:,3) + x(:,5))) + M3.*(L3.*cos(x(:,1) + x(:,5)) - MZ3.*cos(x(:,1) + x(:,5)) - MY3.*sin(x(:,1) + x(:,5)) + L4.*cos(x(:,1) + x(:,3) + x(:,5))) + M3.*(L3.*cos(x(:,1) + x(:,5)) - MZ3.*cos(x(:,2) + x(:,5)) - MY3.*sin(x(:,2) + x(:,5)) + L4.*cos(x(:,1) + x(:,3) + x(:,5))))./(M1 + 2.*M3 + 2.*M4) + rVPP.*cos(gam) + L4.*cos(x(:,1) + x(:,3) + x(:,5));
  out.pHead1 = out.pH1 + -L1*sin(x(5));
  out.pHead2 = out.pH2 + L1*cos(x(5));
  out.pFem11 = out.pH1 + MY3*cos(x(1) + x(5)) - MZ3*sin(x(1) + x(5));
  out.pFem12 = out.pH2 + MZ3*cos(x(1) + x(5)) + MY3*sin(x(1) + x(5));
  out.pFem21 = out.pH1 + MY3*cos(x(2) + x(5)) - MZ3*sin(x(2) + x(5));
  out.pFem22 = out.pH2 + MZ3*cos(x(2) + x(5)) + MY3*sin(x(2) + x(5));
  out.pG11 = out.pH1 + -L3*sin(x(1) + x(5));
  out.pG12 = out.pH2 + L3*cos(x(1) + x(5));
  out.pG21 = out.pH1 + -L3*sin(x(2) + x(5));
  out.pG22 = out.pH2 + L3*cos(x(2) + x(5));
  out.pTib11 = out.pH1 + MY4*cos(x(1) + x(3) + x(5)) - L3*sin(x(1) + x(5)) - MZ4*sin(x(1) + x(3) + x(5));
  out.pTib12 = out.pH2 + L3*cos(x(1) + x(5)) + MZ4*cos(x(1) + x(3) + x(5)) + MY4*sin(x(1) + x(3) + x(5));
  out.pTib21 = out.pH1 + MY4*cos(x(2) + x(4) + x(5)) - L3*sin(x(2) + x(5)) - MZ4*sin(x(2) + x(4) + x(5));
  out.pTib22 = out.pH2 + L3*cos(x(2) + x(5)) + MZ4*cos(x(2) + x(4) + x(5)) + MY4*sin(x(2) + x(4) + x(5));
  out.pFoot11 = out.pH1 + - L3*sin(x(1) + x(5)) - L4*sin(x(1) + x(3) + x(5));
  out.pFoot12 = out.pH2 + L3*cos(x(1) + x(5)) + L4*cos(x(1) + x(3) + x(5));
  out.pFoot21 = out.pH1 + - L3*sin(x(2) + x(5)) - L4*sin(x(2) + x(4) + x(5));
  out.pFoot22 = out.pH2 + L3*cos(x(2) + x(5)) + L4*cos(x(2) + x(4) + x(5));
end
