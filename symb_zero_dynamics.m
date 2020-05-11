function symb_zero_dynamics
%SYMB_ZERO_DYNAMICS    Derive the zero dynamics for the kneed biped
%    walker with Franck Plestan's offset masses model using Bezier
%    polynomials.

%Copyright (c) 2003 by Eric Westervelt and Jessy Grizzle.  This code
%may be freely used for noncommercial ends. If use of this code in
%part or in whole results in publication, proper citation must be
%included in that publication.  This code comes with no guarantees or
%support.

%Eric R. Westervelt & Jessy W. Grizzle
%1/8/01
%11/30/03 - updated for the CDC

if isempty(dir('pm_mat_files_sym/work_symb_zero_dynamics.mat'))
  load pm_mat_files_sym/work_symb_control_5links
  load pm_mat_files_sym/work_symb_model_5links_5DoF
  
  % First step is to do a partial feedback linearization of the
  % dynamics. This will simplify some of the computations and,
  % hopefully, improve the numerical conditioning. See Reyhanoglu,
  % McClamroch and van der Schaft, IEEE T-AC, Sept. 1999, pp. 1663.

  %
  % Change coordinates
  %
  [qr,dqr,q_b,dq_b] = rel2abs;
  %
  q31 = q_b(1);
  q32 = q_b(2);
  q41 = q_b(3);
  q42 = q_b(4);
  q1 =  q_b(5);
  dq31 = dq_b(1);
  dq32 = dq_b(2);
  dq41 = dq_b(3);
  dq42 = dq_b(4);
  dq1 =  dq_b(5);
  q =  qr;
  dq = dqr;
  
  [n,m] = size(B);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Change of coordinates -- attempting to complete with a result of
  % Frobenius's Theorem
  %
  %xh = [H; LfH; theta_s_];
  
  % compute augmented B matrix
  B_aug = [B,[zeros(4,1);1]];
  
  % compute last entry of change of coordinates as last row of ...
  tmp_abs = inv([B [zeros(4,1);1]])*D*dq;
  gamma_s = tmp_abs(5);

  disp('[Computed change of coordinates]');
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Invert map
  %
  syms dz1 dz2 z1 z2 real
  
  syms hd1 hd2 hd3 hd4 real % other terms are zero
  hd_tmp = [hd1; hd2; hd3; hd4];
  
  q_s_zd = phi_inv*[hd_tmp; z1];
  
  q_s_zd = simplify(q_s_zd);
  
  disp('[computed position coordinates restricted to zero dynamics]');
  disp('[manifold                                                 ]');
  
  syms jac_hd1 jac_hd2 jac_hd3 jac_hd4 real
  jac_hd_tmp = [jac_hd1; jac_hd2; jac_hd3; jac_hd4];
  temp_1 = [-jac_hd_tmp*jacobian(theta_s_,q); D(5,:)];
  
  disp('[computed parts needed to compute velocity coordinates]');
  disp('[restricted to zero dynamics manifold                 ]');
  
  dq_s_zd_almost = inv([phi(1:4,:);zeros(1,5)] + temp_1)* ...
      [zeros(4,1); z2];
  
  dq_s_zd = subs(dq_s_zd_almost, q(1), q_s_zd(1));
  dq_s_zd = subs(dq_s_zd, q(2), q_s_zd(2));
  dq_s_zd = subs(dq_s_zd, q(3), q_s_zd(3));
  dq_s_zd = subs(dq_s_zd, q(4), q_s_zd(4));
  dq_s_zd = subs(dq_s_zd, q(5), q_s_zd(5));

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Take time derivative in new coordinates...the system with the
  % change of coordinates and decoupling feedback (minus the last
  % two states)...we know this is the form because this is what
  % feedback linearization does!...
  %
  %syms v1 v2 v3 v4 real % inputs
  %dxh = [xh_5; xh_6; xh_7; xh_8; v1; v2; v3; v4];
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Now for the zero dynamics -- assume all outputs are zero
  %
  disp('[Calculating zero dynamics...]');
  
  dz1 = simplify(jacobian(theta_s_,q)*dq_s_zd)
  dz2 = -G(n);
  syms hip_drag_force
  hip_drag = jacobian(pH(2),q)*hip_drag_force;
  dz2_hip_drag = -G(n)+hip_drag(n);

  dz2 = subs(dz2,q(1),q_s_zd(1));
  dz2 = subs(dz2,q(2),q_s_zd(2));
  dz2 = subs(dz2,q(3),q_s_zd(3));
  dz2 = subs(dz2,q(4),q_s_zd(4));
  dz2 = subs(dz2,q(5),q_s_zd(5));
  dz2 = simplify(dz2)
  dz2_hip_drag = subs(dz2_hip_drag,q(1),q_s_zd(1));
  dz2_hip_drag = subs(dz2_hip_drag,q(2),q_s_zd(2));
  dz2_hip_drag = subs(dz2_hip_drag,q(3),q_s_zd(3));
  dz2_hip_drag = subs(dz2_hip_drag,q(4),q_s_zd(4));
  dz2_hip_drag = subs(dz2_hip_drag,q(5),q_s_zd(5));
  dz2_hip_drag = simplify(dz2_hip_drag);
  disp('[The zero dynamics in z1, z2 coordinates !]');
  
  jac_hd = jacobian(hd,z1);
  jac2_hd = jacobian(jac_hd,z1);
  
  syms jac2_hd1 jac2_hd2 jac2_hd3 jac2_hd4 real
  jac2_hd_tmp = [jac2_hd1; jac2_hd2; jac2_hd3; jac2_hd4];

  alpha = z2/dz1;
  beta = dz2;
  jac_alpha = jacobian(alpha,[hd_tmp;jac_hd_tmp]) ...
      *[jac_hd_tmp; jac2_hd_tmp];
  
  save pm_mat_files_sym/work_symb_zero_dynamics
  disp('[save mat_files_sym/work_symb_zero_dynamics]');
else
  load pm_mat_files_sym/work_symb_zero_dynamics
  disp('[DONE loading mat_files_sym/work_symb_zero_dynamics.mat]')
end
%
%
%
fid=print_function_header('zd_project',...
			  'z',...
			  'x',...
			  ['Project to the zero dynamics manifold.']);

fprintf(fid,'modelP;\n\n');

% Reassign configuration 
fprintf(fid,'[nn,mm] = size(x);\n');
fprintf(fid,'if (nn > 1 & mm == 10)\n');
fprintf(fid,'  q31L=x(:,1); q32L=x(:,2); q41L=x(:,3); q42L=x(:,4); q1L=x(:,5);\n');
fprintf(fid,'  dq31L=x(:,6); dq32L=x(:,7); dq41L=x(:,8); dq42L=x(:,9); dq1L=x(:,10);\n');
fprintf(fid,'else\n');
fprintf(fid,'  q31L=x(1); q32L=x(2); q41L=x(3); q42L=x(4); q1L=x(5);\n');
fprintf(fid,'  dq31L=x(6); dq32L=x(7); dq41L=x(8); dq42L=x(9); dq1L=x(10);\n');
fprintf(fid,'end\n\n');


fprintf(fid,'%% theta (in absolute coordinates)\n');
ttt = vectorize(char(theta_s_));
fprintf(fid,'z1 = %s;\n\n',ttt);

fprintf(fid,'%% gamma (in absolute coordinates)\n');
ttt = vectorize(char(gamma_s));
fprintf(fid,'z2 = %s;\n\n',fixlength(ttt,'*+-',65,'       '));
fprintf(fid,'z3 = 0;\n');

fprintf(fid,'z = [z1 z2 z3];\n');
fclose(fid);
%
%
%
fid=print_function_header('zd',...
			  'dz',...
			  'z',...
			  ['']);

fprintf(fid,'[nn,mm] = size(z);\n');
fprintf(fid,'if (nn > 1 & mm == 2)\n');
fprintf(fid,'  z1 = z(:,1); z2 = z(:,2);\n');
fprintf(fid,'else\n');
fprintf(fid,'  z1 = z(1); z2 = z(2);\n');
fprintf(fid,'end\n\n');

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

% State one
ttt = char(vectorize(dz1));
fprintf(fid,'%% dz1\n');
fprintf(fid,'dz1 = %s;\n\n',fixlength(ttt,'*+-',65,'       '));

% State two
ttt = char(vectorize(dz2));
fprintf(fid,'%% dz2\n');
fprintf(fid,'dz2 = %s;\n\n',fixlength(ttt,'*+-',65,'       '));
fprintf(fid,'x_proj = zd_s_lift(z);\n');
fprintf(fid,'qdot_proj =x_proj(6:9);\n');
fprintf(fid,'u_proj = control_synth(x_proj);\n');
fprintf(fid,'dz3 = sum(u_proj.^2);\n');
fprintf(fid,'dz = [dz1 dz2 dz3];\n');
fclose(fid);
