function symb_model_5links_7DOF
%SYMB_MODEL_LAP

%Copyright (c) 2003 by Eric Westervelt and Jessy Grizzle.  This code
%may be freely used for noncommercial ends. If use of this code in
%part or in whole results in publication, proper citation must be
%included in that publication.  This code comes with no guarantees or
%support.

%Eric R. Westervelt & Jessy W. Grizzle
% Minh edit 2016 

if isempty(dir('pm_mat_files_sym/work_symb_model_5links_7DoF.mat'))
  % This is for a robot with an upright trunk, two legs and knees. 
  % The model is for five degrees of freedom of the robot, 
  % with Foot1 touching the ground, 
  % q1: absolute angle of torso
  % q31: torso-femur 1; q32: torso - femur 2
  % q41: femur 1 - tibia 1; q42: femur 2 - tibia 2
  syms q31L q32L q41L q42L q1L real % relative angle
  syms dq31L dq32L dq41L dq42L dq1L real
  syms y z dy dz ycm zcm dycm dzcm real
  syms ddq31L ddq32L ddq41L ddq42L ddq1L
  syms g L1 L3 L4 M1 M3 M4 real
  syms MY1 MY3 MY4 MZ1 MZ3 MZ4 XX1 XX3 XX4 real
  % Change coordinates % qxx absolute angle; q: relative angle
  [q_rel,dq_rel,q_abs,dq_abs] = rel2abs;
  q31 = q_abs(1);  
  q32 = q_abs(2);
  q41 = q_abs(3);
  q42 = q_abs(4);
  q1 =  q_abs(5);
  dq31 = dq_abs(1);
  dq32 = dq_abs(2);
  dq41 = dq_abs(3);
  dq42 = dq_abs(4);
  dq1 =  dq_abs(5); 
  q =  q_rel;
  dq = dq_rel;
  % add ycm and zcm
  q  = [q; ycm; zcm];
  q_old = [q(1:5); y; z];
  dq = [dq; dycm; dzcm];
  dq_old = [dq(1:5); dy; dz];
  %
  % L3: Length of femur; L4: Length of Tibia
  R = @(t) [cos(t) -sin(t); sin(t) cos(t)]; % rotation 
  pH    = [y; z];
  pT    = pH + R(q1)*[MY1; MZ1];
  pFem1 = pH + R(q31)*[MY3; MZ3];
  pFem2 = pH + R(q32)*[MY3; MZ3];
  pG1   = pH + R(q31)*[0;1]*L3;
  pG2   = pH + R(q32)*[0;1]*L3;
  pTib1 = pG1 + R(q41)*[MY4; MZ4];
  pTib2 = pG2 + R(q42)*[MY4; MZ4];
  %
  %
  M_total = M1 + 2*M3 + 2*M4;
  cm = simplify(1/M_total*(M1*pT+M3*pFem1+M3*pFem2+M4*pTib1+M4*pTib2));
  dcm = simplify(jacobian(cm,q_old)*dq_old);
  %
  y = solve(cm(1)-ycm,'y'); y = simplify(y); 
  z = solve(cm(2)-zcm,'z'); z = simplify(z);
  dy = solve(dcm(1)-dycm,'dy'); dy = simplify(dy);
  dz = solve(dcm(2)-dzcm,'dz'); dz = simplify(dz);
  %
  pH = [y; z];
  pT    = pH + R(q1)*[MY1; MZ1];
  pFem1 = pH + R(q31)*[MY3; MZ3];
  pFem2 = pH + R(q32)*[MY3; MZ3];
  pG1   = pH + R(q31)*[0;1]*L3;
  pG2   = pH + R(q32)*[0;1]*L3;
  pTib1 = pG1 + R(q41)*[MY4; MZ4];
  pTib2 = pG2 + R(q42)*[MY4; MZ4];
  pFoot1= pG1 + R(q41)*[0;1]*L4;
  pFoot2= pG2 + R(q42)*[0;1]*L4;
  pHead = pH +  R(q1)*[0;1]*L1;
  % velocities
  vT=jacobian(pT,q)*dq;
  vFem1 = jacobian(pFem1,q)*dq;
  vFem2 = jacobian(pFem2,q)*dq;
  vTib1 = jacobian(pTib1,q)*dq;
  vTib2 = jacobian(pTib2,q)*dq;
  vFoot1=jacobian(pFoot1,q)*dq;
  vFoot2=jacobian(pFoot2,q)*dq;
  %
  KET    = simplify(1/2*M1*vT.'*vT + 1/2*XX1*dq1^2);
  % actuator interia only depends on the relative link angle (see
  % the IAx term)
  KEFem1 = simplify(1/2*M3*vFem1.'*vFem1 + 1/2*XX3*dq31^2);
  KEFem2 = simplify(1/2*M3*vFem2.'*vFem2 + 1/2*XX3*dq32^2);
  KETib1 = simplify(1/2*M4*vTib1.'*vTib1 + 1/2*XX4*dq41^2);
  KETib2 = simplify(1/2*M4*vTib2.'*vTib2 + 1/2*XX4*dq42^2);
		   
  KE = simplify(KET + KEFem1 + KEFem2 + KETib1 + KETib2);
  %
  PE = g*(M1*pT(2) + M3*pFem1(2) + M3*pFem2(2) + M4*pTib1(2) + M4*pTib2(2));
	  
  PE = simplify(PE);
  %
  %
  % Model NOTATION: Spong and Vidyasagar, page 142, Eq. (6.3.12)
  %                 D(q)ddq + C(q,dq)*dq + G(q) = B*tau + E1*F_external
  %
  %L=KE-PE;
  G = jacobian(PE,q).';
  G = simplify(G);
  D = jacobian(KE,dq).';
  n = max(size(q));
  for k = 1:n
    D_tmp = D(k);
    tic
    D_tmp2 = simplify(D_tmp);
    toc
    D(k) = D_tmp2;
    disp(['D simplify ',num2str(k)]);
  end
  D = jacobian(D,dq);
  for k = 1:n
    for j = 1:k
      D(k,j) = simplify(D(k,j));
      disp(['D simplify (',num2str(k),',',num2str(j),')']);
    end
    D(j,k) = D(k,j);
  end
  %
  C = sym(zeros(length(q), length(q)));
  n = max(size(q));
  for k = 1:n
    for j = 1:n
      C(k,j) = 0*g;
      for i=1:n
      C(k,j) = C(k,j)+1/2*(diff(D(k,j),q(i)) + ...
			     diff(D(k,i),q(j)) - ...
			     diff(D(i,j),q(k)))*dq(i);
      end
    end
  end
  C = simplify(C);
  % C = simplify(C);
  %
  % Compute the matrix for the input torques and then the contact
  % points.
  %
  Phi_0 = [q31-q1; q32-q1; q41-q31; q42-q32; q1; y; z];
  B = jacobian(Phi_0,q);
  B = B.'*[eye(4,4); zeros(3,4)];
  %
  %
  % F_ext = [F_T^1;F_N^1;F_T^2;F_N^2];
  %
  Phi_1 = [pFoot1; pFoot2];
  E1 = jacobian(Phi_1,q).';
  %
  % Compute matrices associated with the viscous and static
  % friction
  %
  %
  % Now, compute the forces acting on the support foot, assuming that
  % the robot is in the swing phase
  %
  c = pFoot1;
  jac_c = jacobian(c,q);             % c_dot  = jac_c * dq
  lfc = jac_c*dq;
  g1 = jac_c.';
  gamma_c = jacobian(jac_c*dq,q)*dq; % c_ddot = gamma_c + jac_c * ddq
  
  save pm_mat_files_sym/work_symb_model_5links_7DoF
  disp('[DONE saving MAT_FILES_SYM/work_symb_model_5links_7DoF.MAT]')
else
  load pm_mat_files_sym/work_symb_model_5links_7DoF
  disp('[DONE loading MAT_FILES_SYM/work_symb_model_5links_7DoF.MAT]')
end
N = max(size(q));
fid = fopen('m_files/modelP.m','w');
print_modelParameters(fid);
fclose(fid);

%%
fid=print_function_header('dyn_mod_5links_7DoF',...
			  'D,C,G,B,E1',...
			  'q,dq',...
			  ['is the kneed biped walking model.']);
fprintf(fid,'modelP;\n\n');
fprintf(fid,['\nq31L=q(1);q32L=q(2);q41L=q(3);q42L=q(4);q1L=q(5);' ...
	       'ycm=q(6);zcm=q(7);']);
fprintf(fid,['\ndq31L=dq(1);dq32L=dq(2);dq41L=dq(3);dq42L=dq(4);' ...
	       'dq1L=dq(5);dycm=dq(6);dzcm=dq(7);']);
%%%%%
fprintf(fid,'%% D matrix\n');
fprintf(fid,'D=zeros(%s);\n',num2str(N));
for k=1:N,
  for j=1:N,
    if D(k,j)~=0
      ttt=fixlength(char(D(k,j)),'*+-',65,'       ');
      fprintf(fid,'D(%s,%s)=%s;\n',num2str(k),num2str(j),ttt);
    end
  end
end

fprintf(fid,'\n%% C matrix\n');
fprintf(fid,'C=zeros(%s);\n',num2str(N));
for k=1:N,
  for j=1:N,
    if C(k,j)~=0
      ttt=fixlength(char(C(k,j)),'*+-',65,'       ');    
      fprintf(fid,'C(%s,%s)=%s;\n',num2str(k),num2str(j),ttt);
    end
  end
end

fprintf(fid,'\n%% G matrix\n');
fprintf(fid,'G=zeros(%s,1);\n',num2str(N));
for k=1:N,
  if G(k)~=0
    ttt=fixlength(char(G(k)),'*+-',65,'     ');
    fprintf(fid,'G(%s)=%s;\n',num2str(k),ttt);
  end
end

fprintf(fid,'\n%% B matrix\n');
[N,M]=size(B);
fprintf(fid,'B=zeros(%s,%s);\n',num2str(N),num2str(M));
for k=1:N,
  for j=1:M,
    if B(k,j)~=0
      ttt=char(B(k,j));
      fprintf(fid,'B(%s,%s)=%s;\n',num2str(k),num2str(j),ttt);
    end
  end
end
[n,m]=size(E1);
fprintf(fid,'\n%s',['E1=zeros(',num2str(n),',',num2str(m),');']);
for i=1:n
  for j=1:m
    Temp0=E1(i,j);
    if Temp0~=0
      Temp1=char(Temp0);
      Temp2=['E1(',num2str(i),',',num2str(j),')=',Temp1,';'];
      fprintf(fid,'\n%s',Temp2);
    end
  end
end
fclose(fid);

%%
%print_model
fid=print_function_header('stance_foot_vel',...
			  'vh,vv',...
			  'x',...
			  ['']);
fprintf(fid,'modelP;\n\n');
fprintf(fid,['q31L=x(1);q32L=x(2);q41L=x(3);q42L=x(4);' ...
	     'q1L=x(5); ycm=x(6);zcm=x(7);\n']);
fprintf(fid,['dq31L=x(8);dq32L=x(9);dq41L=x(10);dq42L=x(11);' ...
	     'dq1L=x(12); dycm=x(13);dzcm=x(14);\n\n']);
ttt = char(vFoot1(1));
fprintf(fid,'vh = %s;\n',fixlength(ttt,'*+-',65,'    '));
ttt = char(vFoot1(2));
fprintf(fid,'vv = %s;\n',fixlength(ttt,'*+-',65,'    '));
fclose(fid);
