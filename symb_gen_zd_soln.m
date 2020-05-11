function symb_gen_zd_soln
%SYMB_GEN_ZD_SOLN   Derive the zero dynamics solution.

%Copyright (c) 2003 by Eric Westervelt and Jessy Grizzle.  This code
%may be freely used for noncommercial ends. If use of this code in
%part or in whole results in publication, proper citation must be
%included in that publication.  This code comes with no guarantees or
%support.

%Eric R. Westervelt & Jessy W. Grizzle
%8/6/01
%11/30/03 - updated for the CDC

if isempty(dir('pm_mat_files_sym/work_symb_zd_soln.mat'))
  if isempty(dir('pm_mat_files_sym/work_symb_zero_dynamics.mat'))
    error('run symb_zero_dynamics.m')
  else
    load pm_mat_files_sym/work_symb_zero_dynamics
  end
  syms z1
  alpha = z2/dz1;
  beta = dz2;
  inv_dz1 =1/dz1;
  
  save pm_mat_files_sym/work_symb_zd_soln
else
  load pm_mat_files_sym/work_symb_zd_soln
end
%
%
%
fid=print_function_header('alpha',...
			  'alpha',...
			  'z1,a,b,m',...
			  ['']);
fprintf(fid,'modelP;\n\n');
print_controlParameters(fid,M);

for k = 1:4
  ttt = char(vectorize(hd(k)));
  fprintf(fid,'hd%d = %s;\n',k,fixlength(ttt,'*+-',65,'      '));
end
fprintf(fid,'\n');
for k = 1:4
  ttt = char(vectorize(jac_hd(k)));
  fprintf(fid,'jac_hd%d = %s;\n',k,fixlength(ttt,'*+-',65,'          '));
end
fprintf(fid,'\n');

ttt = fixlength(vectorize(char(alpha)),'*+-',65,'        ');
fprintf(fid,'alpha = %s;',ttt);
fclose(fid);
%
%
%
fid=print_function_header('beta',...
			  'beta',...
			  'z1,a,b,m',...
			  ['']);
fprintf(fid,'modelP;\n\n');
print_controlParameters(fid,M);

for k = 1:4
  ttt = char(vectorize(hd(k)));
  fprintf(fid,'hd%d = %s;\n',k,fixlength(ttt,'*+-',65,'      '));
end
fprintf(fid,'\n');
for k = 1:4
  ttt = char(vectorize(jac_hd(k)));
  fprintf(fid,'jac_hd%d = %s;\n',k,fixlength(ttt,'*+-',65,'          '));
end
fprintf(fid,'\n');

ttt = fixlength(vectorize(char(beta)),'*+-',65,'       ');
fprintf(fid,'beta = %s;',ttt);
fclose(fid);
%
%
%
fid=print_function_header('alpha_beta',...
			  'alpha_times_beta',...
			  'z1,a,b,m',...
			  ['']);
fprintf(fid,'alpha_val = alpha(z1,a,b,m);\n');
fprintf(fid,'beta_val = beta(z1,a,b,m);\n\n');
fprintf(fid,'alpha_times_beta = alpha_val.*beta_val;');
fclose(fid);
%
%
fid=print_function_header('inv_dz1',...
			  'inv_dz1',...
			  'z1,a,b,m',...
			  ['']);
fprintf(fid,'modelP;\n\n');
print_controlParameters(fid,M);

for k = 1:4
  ttt = char(vectorize(hd(k)));
  fprintf(fid,'hd%d = %s;\n',k,fixlength(ttt,'*+-',65,'      '));
end
fprintf(fid,'\n');
for k = 1:4
  ttt = char(vectorize(jac_hd(k)));
  fprintf(fid,'jac_hd%d = %s;\n',k,fixlength(ttt,'*+-',65,'          '));
end
fprintf(fid,'\n');

ttt = fixlength(vectorize(char(inv_dz1)),'*+-',65,'       ');
fprintf(fid,'inv_dz1 = %s;',ttt);
fclose(fid);
%
%
%