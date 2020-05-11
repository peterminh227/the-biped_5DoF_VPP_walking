function print_controlParameters_jacH(fid,M)
%PRINT_CONTROLPARAMETERS

%Copyright (c) 2003 by Eric Westervelt and Jessy Grizzle.  This code
%may be freely used for noncommercial ends. If use of this code in
%part or in whole results in publication, proper citation must be
%included in that publication.  This code comes with no guarantees or
%support.

%Eric R. Westervelt & Jessy W. Grizzle
%11/30/03 - updated for the CDC
% 2016, P.M modified for optimation
%fprintf(fid,'[a,b,m] = controlParameters;\n\n');

a = char('a');
for j = 0:3
  for k = 0:M
    fprintf(fid,[char(a+j),num2str(k),'=a(', ...
		 num2str(1+k+j*(M+1)),'); ']);
  end
  fprintf(fid,'\n');
end
fprintf(fid,'\n');

