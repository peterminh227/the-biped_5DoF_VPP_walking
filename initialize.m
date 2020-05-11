%INITIALIZE.M
%
%   Initialize the biped robot

%Copyright (c) 2003 by Eric Westervelt and Jessy Grizzle.  This code
%may be freely used for noncommercial ends. If use of this code in
%part or in whole results in publication, proper citation must be
%included in that publication.  This code comes with no guarantees or
%support.

%Eric R. Westervelt & Jessy W. Grizzle
%1/6/01
%11/30/03 - updated for the CDC
% > Peter Minh fix
%[a,b,m,dtheta_fixed] = controlParameters;
function x0 = initialize(optim_var)
x0 = sigma(optim_var);
out = impact_5links_5DoF(x0);
% return 7 DoF
%out = red2full(out)';
x0 = out(1:10);
end
