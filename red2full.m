function x_com = red2full(x_red)
%RED2FULL

%Copyright (c) 2003 by Eric Westervelt and Jessy Grizzle.  This code
%may be freely used for noncommercial ends. If use of this code in
%part or in whole results in publication, proper citation must be
%included in that publication.  This code comes with no guarantees or
%support.

%Eric R. Westervelt & Jessy W. Grizzle
%5/7/02
%11/30/03 - updated for the CDC
% P.M update 2016
[n,m] = size(x_red);

[cmh,cmv,dcmh,dcmv] = center_of_mass_5DoF(x_red);

if n ~= 1 && m ~= 1
  x_com = [x_red(:,1:5) cmh cmv x_red(:,6:10) dcmh dcmv];
else
  if n < m
    x_com = [x_red(1:5) cmh cmv x_red(6:10) dcmh dcmv];
  else
    x_com = [x_red(1:5); cmh; cmv; x_red(6:10); dcmh; dcmv];
  end
end