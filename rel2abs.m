function [q_rel,dq_rel,q_abs,dq_abs] = rel2abs(x_abs)
%rel2abs   Change from my absolute coorinates to the Lapin
%    coordinates.

%Copyright (c) 2003 by Eric Westervelt and Jessy Grizzle.  This code
%may be freely used for noncommercial ends. If use of this code in
%part or in whole results in publication, proper citation must be
%included in that publication.  This code comes with no guarantees or
%support.

if nargin == 0
  syms q31L q32L q41L q42L q1L real
  syms dq31L dq32L dq41L dq42L dq1L real
  q_rel  = [q31L;  q32L;  q41L;  q42L;  q1L];
  dq_rel = [dq31L; dq32L; dq41L; dq42L; dq1L];
end

Trel2abs = [ 1   0   0   0   1;
             0   1   0   0   1;
             1   0   1   0   1;
             0   1   0   1   1;
             0   0   0   0   1];

brel2abs = [0; 0; 0; 0; 0];

%q_abs  = Trel2abs*q_lap + brel2abs;
%dq_abs = Trel2abs*dq_lap;
%
%q_rel  = inv(Trel2abs)*q_abs - brel2abs;
%dq_rel = inv(Trel2abs)*dq_abs;
if nargin == 0
  % relative coordinates in terms of actual coords
  q_abs  = Trel2abs*q_rel + brel2abs; 
  dq_abs = Trel2abs*dq_rel;
else
  q_rel(1:5)  = inv(Trel2abs)*(x_abs(1:5) - brel2abs);
  q_rel(6:10) = Trel2abs*x_abs(6:10);
  q_rel = q_rel';
  q_rel_sym = 0;
  dq_rel_sym = 0;
end
