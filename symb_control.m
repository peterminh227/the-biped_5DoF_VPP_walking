function symb_control
%SYMB_CONTROL.M

%Copyright (c) 2003 by Eric Westervelt and Jessy Grizzle.  This code
%may be freely used for noncommercial ends. If use of this code in
%part or in whole results in publication, proper citation must be
%included in that publication.  This code comes with no guarantees or
%support.

%Eric R. Westervelt & Jessy W. Grizzle
%2/15/01
%11/30/03 - updated for the CDC

if isempty(dir('pm_mat_files_sym/work_symb_control_5links.mat'))
  load pm_mat_files_sym/work_symb_model_5links_5DoF

  syms u z1 dz1 ddz1 real
  
  syms a0 a1 a2 a3 a4 a5 a6
  syms b0 b1 b2 b3 b4 b5 b6
  syms c0 c1 c2 c3 c4 c5 c6
  syms d0 d1 d2 d3 d4 d5 d6
  
  a = [a0 a1 a2 a3 a4 a5 a6];
  b = [b0 b1 b2 b3 b4 b5 b6];
  c = [c0 c1 c2 c3 c4 c5 c6];
  d = [d0 d1 d2 d3 d4 d5 d6];
  
  hd1 = 0; hd2 = 0; hd3 = 0; hd4 = 0;
  
  % form bezier parameterization of polynomials
  M = 6;
  for k = 1:M+1
    kk = k - 1;
    hd1 = hd1 + a(k)*factorial(M)/ ...
	 factorial(kk)/factorial(M-kk)*u^(kk)*(1-u)^(M-kk);
    hd2 = hd2 + b(k)*factorial(M)/ ...
	 factorial(kk)/factorial(M-kk)*u^(kk)*(1-u)^(M-kk);
    hd3 = hd3 + c(k)*factorial(M)/ ...
	 factorial(kk)/factorial(M-kk)*u^(kk)*(1-u)^(M-kk);
    hd4 = hd4 + d(k)*factorial(M)/ ...
	 factorial(kk)/factorial(M-kk)*u^(kk)*(1-u)^(M-kk);
  end
  
  syms b m
  
  hd1 = subs(hd1,u,(z1-b)/m);
  hd2 = subs(hd2,u,(z1-b)/m);
  hd3 = subs(hd3,u,(z1-b)/m);
  hd4 = subs(hd4,u,(z1-b)/m);
  
  H0 = [1  0  0  0  0;
		0  1  0  0  0;
		0  0  1  0  0;
		0  0  0  1  0];
  
  hd = [hd1;
		hd2;
		hd3;
		hd4];
  
  % The output function:
  %
  %   a linear term, constant term, plus a nonlinear term
  %   parameterized by theta
  %
  h = H0*q-hd;
  theta_s_  = -q41L/2 - q31L - q1L;
  dtheta_s_ = -dq41L/2 - dq31L - dq1L;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Invert position coordinates
  %
  syms hd1 hd2 hd3 hd4 real
  syms jac_hd1 jac_hd2 jac_hd3 jac_hd4 real
  syms z1 z2 real
  
  hd_tmp = [hd1; hd2; hd3; hd4];
  jac_hdtmp = [jac_hd1; jac_hd2; jac_hd3; jac_hd4];
  jac_hd  = jacobian(hd,z1);
  jac2_hd = jacobian(jac_hd,z1);
  
  phi = [H0;
	 jacobian(theta_s_,q)];
 
  phi_inv = inv(phi);
  
  q_rel_zd = phi_inv*[hd_tmp; z1];
  q_rel_zd = simplify(q_rel_zd);

  % this is gamma_f0
  D_5 = subs(D(5,:), {q31L, q32L, q41L, q42L, q1L}, ...
	     {q_rel_zd(1), q_rel_zd(2), q_rel_zd(3), q_rel_zd(4), q_rel_zd(5)});
  temp1 = phi(1:4,:) - jac_hdtmp*jacobian(theta_s_,q_rel);
  temp2 = inv([temp1;
	       D_5]);

  dq_rel_zd = temp2*[zeros(4,1); z2];

  
  jac_h = H0 - jacobian(hd,z1)*jacobian(theta_s_,q);

  lfh = H0*dq - jacobian(hd,z1)*dz1;
  
  % y_ddot = gamma + jac_h * ddq
  gamma = jacobian(jacobian(-hd,z1),z1)*dz1^2; 
  
  obj_fun_foot_on_gnd = subs(pFoot2(2), ...
			     {q31L, q32L, q41L, q42L, q1L}, ...
			     {q_rel_zd(1), q_rel_zd(2), q_rel_zd(3), ...
		    q_rel_zd(4),q_rel_zd(5)});

  hd_end = subs(hd,z1,b+m);
  jac_hd_end = subs(jac_hd,z1,b+m);
  
  obj_fun_foot_on_gnd = subs(obj_fun_foot_on_gnd, ...
			     {hd1, hd2, hd3, hd4}, ...
			     {hd_end(1), hd_end(2), hd_end(3), ...
		    hd_end(4)});
  
  sigma_pos = subs(q_rel_zd, ...
		   {hd1, hd2, hd3, hd4, z1}, ...
		   {hd_end(1), hd_end(2), hd_end(3), hd_end(4), b+m});
  sigma_vel = subs(dq_rel_zd, ...
		   {hd1, hd2, hd3, hd4, ...
		    jac_hd1, jac_hd2, jac_hd3, jac_hd4, z1}, ...
		   {hd_end(1), hd_end(2), hd_end(3), hd_end(4), ...
		    jac_hd_end(1), jac_hd_end(2), jac_hd_end(3), ...
		    jac_hd_end(4), b+m});
  
  save pm_mat_files_sym/work_symb_control_5links
  disp('[WROTE work_symb_control_5links.mat]')
else
  load pm_mat_files_sym/work_symb_control_5links
  disp('[DONE loading mat_files_sym/work_symb_control_5links.mat]')
end


%------------------------------------------------------
%
% Automatically generate .m-files needed for simulation
%
%> modify for optimization zero dynamics. Peter Minh 2016
fid=print_function_header('decouple_5links_5DoF',...
			  'h,lfh,l2fh,lglfh,jac_h,inv_D,C,G,B',...
			  'q,dq,a,b,m',...
			  ['']);

fprintf(fid,'q31L=q(1);q32L=q(2);q41L=q(3);q42L=q(4);q1L=q(5);');
fprintf(fid,'\ndq31L=dq(1);dq32L=dq(2);dq41L=dq(3);');
fprintf(fid,'dq42L=dq(4);dq1L=dq(5);\n\n');

print_controlParameters(fid,M);

fprintf(fid,'\ns = (%s - b)/m;\n',char(theta_s_));
fprintf(fid,'if s > 1\n');
fprintf(fid,'  s = 1;\n');
fprintf(fid,'elseif s < 0\n');
fprintf(fid,'  s = 0;\n');
fprintf(fid,'end\n');
fprintf(fid,'z1 = m*s + b;\n\n');

fprintf(fid,'dz1 = %s;\n\n',char(dtheta_s_));

p = max(size(h));
fprintf(fid,['h = zeros(',num2str(p),',1);']);
for i = 1:p
  Temp = char(h(i));
  fprintf(fid,'\nh(%d) = %s;', ...
	  i,fixlength(Temp,'*+-',63,'       '));
end
fprintf(fid,'\n%%\n');
p=max(size(lfh));
fprintf(fid,['lfh = zeros(',num2str(p),',1);\n']);
for i = 1:p
  Temp = char(lfh(i));
  fprintf(fid,'lfh(%d) = %s;\n', ...
	  i,fixlength(Temp,'*+-',63,'         '));
end
fprintf(fid,'\n%%\n');
[p,n] = size(jac_h);
fprintf(fid,['jac_h = zeros(',num2str(p),',',num2str(n),');\n']);
for i = 1:p
  for j = 1:n
    Temp0 = jac_h(i,j);
    if Temp0 ~= 0
      Temp = char(Temp0);
      fprintf(fid,'\njac_h(%d,%d) = %s;', ...
	      i,j,fixlength(Temp,'*+-',60,'             '));
    end
  end
end
fprintf(fid,'\n%%');
fprintf(fid,'\n%%');
fprintf(fid,'\n%%');
p = max(size(gamma));
fprintf(fid,['\ngamma = zeros(',num2str(p),',1);']);
for i = 1:p
  Temp = char(gamma(i));
  fprintf(fid,'\ngamma(%d) = %s;', ...
	  i,fixlength(Temp,'*+-',63,'           '));
end
fprintf(fid,'\n%%');
fprintf(fid,'\n%%');
fprintf(fid,'\n[D,C,G,B] = dyn_mod_5links_5DoF(q,dq);');
fprintf(fid,'\ninv_D = inv(D);');
fprintf(fid,'\nl2fh  = gamma-jac_h*inv_D*(C*dq+G);');
fprintf(fid,'\nlglfh = jac_h*inv_D*B;');
fclose(fid);
%
%
fid=print_function_header('jacobianH',...
			  'jac_h',...
			  'q,dq,a,b,m',...
			  ['']);

fprintf(fid,'q31L=q(1);q32L=q(2);q41L=q(3);q42L=q(4);q1L=q(5);');
fprintf(fid,'\ndq31L=dq(1);dq32L=dq(2);dq41L=dq(3);');
fprintf(fid,'dq42L=dq(4);dq1L=dq(5);\n\n');

print_controlParameters_jacH(fid,M);

fprintf(fid,'\ns = (%s - b)/m;\n',char(theta_s_));
fprintf(fid,'if s > 1\n');
fprintf(fid,'  s = 1;\n');
fprintf(fid,'elseif s < 0\n');
fprintf(fid,'  s = 0;\n');
fprintf(fid,'end\n');
fprintf(fid,'z1 = m*s + b;\n\n');

fprintf(fid,'dz1 = %s;\n\n',char(dtheta_s_));


[p,n] = size(jac_h);
fprintf(fid,['jac_h = zeros(',num2str(p),',',num2str(n),');\n']);
for i = 1:p
  for j = 1:n
    Temp0 = jac_h(i,j);
    if Temp0 ~= 0
      Temp = char(Temp0);
      fprintf(fid,'\njac_h(%d,%d) = %s;', ...
	      i,j,fixlength(Temp,'*+-',60,'             '));
    end
  end
end
fclose(fid);
%
fid=print_function_header('zd_s_lift',...
			  'x',...
			  'z',...
			  ['Lifts the zero dynamics to the full ',...
		    'state of the system.']);

fprintf(fid,'modelP;\n\n');
print_controlParameters(fid,M);

fprintf(fid,'[nn,mm] = size(z);\n');
fprintf(fid,'if (nn > 1 & mm == 2)\n');
fprintf(fid,'  z1 = z(:,1); z2 = z(:,2);\n');
fprintf(fid,'else\n');
fprintf(fid,'  z1 = z(1); z2 = z(2);\n');
fprintf(fid,'end\n\n');
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

ttt = char(vectorize(q_rel_zd(1)));
fprintf(fid,'q31L = %s;\n\n',fixlength(ttt,'*+-',65,'      '));
ttt = char(vectorize(q_rel_zd(2)));
fprintf(fid,'q32L = %s;\n\n',fixlength(ttt,'*+-',65,'      '));
ttt = char(vectorize(q_rel_zd(3)));
fprintf(fid,'q41L = %s;\n\n',fixlength(ttt,'*+-',65,'      '));
ttt = char(vectorize(q_rel_zd(4)));
fprintf(fid,'q42L = %s;\n\n',fixlength(ttt,'*+-',65,'      '));
ttt = char(vectorize(q_rel_zd(5)));
fprintf(fid,'q1L = %s;\n\n',fixlength(ttt,'*+-',65,'     '));
ttt = char(vectorize(dq_rel_zd(1)));
fprintf(fid,'dq31L = %s;\n\n',fixlength(ttt,'*+-',65,'       '));
ttt = char(vectorize(dq_rel_zd(2)));
fprintf(fid,'dq32L = %s;\n\n',fixlength(ttt,'*+-',65,'       '));
ttt = char(vectorize(dq_rel_zd(3)));
fprintf(fid,'dq41L = %s;\n\n',fixlength(ttt,'*+-',65,'       '));
ttt = char(vectorize(dq_rel_zd(4)));
fprintf(fid,'dq42L = %s;\n\n',fixlength(ttt,'*+-',65,'       '));
ttt = char(vectorize(dq_rel_zd(5)));
fprintf(fid,'dq1L = %s;\n\n',fixlength(ttt,'*+-',65,'      '));

fprintf(fid,['x = [q31L q32L q41L q42L q1L dq31L dq32L dq41L dq42L' ...
	     ' dq1L];']);
fclose(fid);
%
%
%
fcn_name = 'sigma';
VAL = b+m;
print_sigma;
%
%
%
%
fid=print_function_header('theta_s',...
			  'z1',...
			  'x',...
			  ['']);

ttt = char(theta_s_);
fprintf(fid,'[nn,mm] = size(x);\n\n');
fprintf(fid,'if nn ~= 1 & mm ~= 1\n');
ttt1 = rep_str(ttt,1);
fprintf(fid,'  z1 = %s;\n',ttt1);
fprintf(fid,'else\n');
ttt2 = rep_str(ttt,2);
fprintf(fid,'  z1 = %s;\n',ttt2);
fprintf(fid,'end\n');
fclose(fid);
%
%
%
fid=print_function_header('dtheta_s',...
			  'dz1',...
			  'x',...
			  ['']);

ttt = char(dtheta_s_);
fprintf(fid,'[nn,mm] = size(x);\n\n');
fprintf(fid,'if nn ~= 1 & mm ~= 1\n');
ttt1 = rep_str(ttt,1);
fprintf(fid,'  dz1 = %s;\n',ttt1);
fprintf(fid,'else\n');
ttt2 = rep_str(ttt,2);
fprintf(fid,'  dz1 = %s;\n',ttt2);
fprintf(fid,'end\n');
fclose(fid);

%------------------------------------------------------

function ttt_out = rep_str(ttt,type)

if type ==1
  ttt = strrep(ttt,'ddq31L','ddq(:,1)');
  ttt = strrep(ttt,'ddq32L','ddq(:,2)');
  ttt = strrep(ttt,'ddq41L','ddq(:,3)');
  ttt = strrep(ttt,'ddq42L','ddq(:,4)');
  ttt = strrep(ttt,'ddq1L', 'ddq(:,5)');
  ttt = strrep(ttt,'dq31L','x(:,6)');
  ttt = strrep(ttt,'dq32L','x(:,7)');
  ttt = strrep(ttt,'dq41L','x(:,8)');
  ttt = strrep(ttt,'dq42L','x(:,9)');
  ttt = strrep(ttt,'dq1L', 'x(:,10)');
  ttt = strrep(ttt,'q31L','x(:,1)');
  ttt = strrep(ttt,'q32L','x(:,2)');
  ttt = strrep(ttt,'q41L','x(:,3)');
  ttt = strrep(ttt,'q42L','x(:,4)');
  ttt = strrep(ttt,'q1L', 'x(:,5)');
  ttt=vectorize(ttt);
else
  ttt = strrep(ttt,'ddq31L','ddq(1)');
  ttt = strrep(ttt,'ddq32L','ddq(2)');
  ttt = strrep(ttt,'ddq41L','ddq(3)');
  ttt = strrep(ttt,'ddq42L','ddq(4)');
  ttt = strrep(ttt,'ddq1L', 'ddq(5)');
  ttt = strrep(ttt,'dq31L','x(6)');
  ttt = strrep(ttt,'dq32L','x(7)');
  ttt = strrep(ttt,'dq41L','x(8)');
  ttt = strrep(ttt,'dq42L','x(9)');
  ttt = strrep(ttt,'dq1L', 'x(10)');
  ttt = strrep(ttt,'q31L','x(1)');
  ttt = strrep(ttt,'q32L','x(2)');
  ttt = strrep(ttt,'q41L','x(3)');
  ttt = strrep(ttt,'q42L','x(4)');
  ttt = strrep(ttt,'q1L', 'x(5)');
end

ttt_out = ttt;
