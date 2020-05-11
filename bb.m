function [out]=bb(u,e_scl)
%BB
%
% u = (y,dy)
%
% out = (nu,v)

%Copyright (c) 2003 by Eric Westervelt and Jessy Grizzle.  This code
%may be freely used for noncommercial ends. If use of this code in
%part or in whole results in publication, proper citation must be
%included in that publication.  This code comes with no guarantees or
%support.

%Eric R. Westervelt & Jessy W. Grizzle
%11/31/03 - updated for the CDC

alpha   = 0.9;
epsilon = 0.075*e_scl;
%
y=u(1);
dy=u(2);
%
dy = epsilon*dy;
phi = y+1/(2-alpha)*sign(dy)*abs(dy)^(2-alpha);
nu = (-sign(dy)*abs(dy)^alpha-sign(phi)*abs(phi)^(alpha/(2-alpha)))/epsilon^2;
%
sigma=1/2;
rho=2;
%
V = (2-alpha)/(3-alpha)*abs(phi)^((3-alpha)/(2-alpha)) + sigma*dy*phi ...
    +rho/(3-alpha)*abs(dy)^(3-alpha);
if V < -1
   nu=0;
end
out=[nu,V];