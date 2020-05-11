%PRINT_SIGMA.M

%Copyright (c) 2003 by Eric Westervelt and Jessy Grizzle.  This code
%may be freely used for noncommercial ends. If use of this code in
%part or in whole results in publication, proper citation must be
%included in that publication.  This code comes with no guarantees or
%support.

%Eric R. Westervelt & Jessy W. Grizzle
%11/30/03 - updated for the CDC

fid=print_function_header(fcn_name,...
			  'x0',...
			  'z2',...
			  ['Maps horizontal velocity of hips ',...
		    'just before impact to state of the system',...
		    ' just before impact.']);
%
fprintf(fid,'modelP;\n\n');
print_controlParameters(fid,M);
%
for k = 1:4
  hd_tmp(k) = subs(hd(k),z1, VAL);
  hd_tmp(k) = simplify(hd_tmp(k));
  fprintf(fid,'hd%d = %s;\n',k,char(hd_tmp(k)));
end
fprintf(fid,'\n');
for k = 1:4
  jac_hd_tmp(k) = subs(jac_hd(k),z1,VAL);
  jac_hd_tmp(k) = simplify(jac_hd_tmp(k));
  fprintf(fid,'jac_hd%d = %s;\n',k,char(jac_hd_tmp(k)));
end
fprintf(fid,'\n');
%
q_rel_zd_tmp = subs(q_rel_zd,z1,VAL);
dq_rel_zd_tmp = subs(dq_rel_zd,z1,VAL);
for k = 1:length(q_rel_zd)
  fprintf(fid,'q(%d) = %s;\n',k, ...
	  fixlength(char(q_rel_zd_tmp(k)),'*+-',65,'        '));
end
for k = 1:length(dq_rel_zd)
  fprintf(fid,'dq(%d) = %s;\n',k, ...
	  fixlength(char(dq_rel_zd_tmp(k)),'*+-',65,'        '));
end
fprintf(fid,'\n');
%
fprintf(fid,'x0 = [q(1:5)''; dq(1:5)''];');
fclose(fid);
