function symb_model_5links_5DOF
% from SYMB_MODEL_LAP_RED
%Copyright (c) 2003 by Eric Westervelt and Jessy Grizzle.  This code
%may be freely used for noncommercial ends. If use of this code in
%part or in whole results in publication, proper citation must be
%included in that publication.  This code comes with no guarantees or
%support.
%From J.W. Grizzle and F. Plestan
%10/27/16 - update for DS phase- P.Minh
if isempty(dir('pm_mat_files_sym/work_symb_model_5links_5DoF.mat'))
  % This is for a robot with an upright trunk, two legs and knees. 
  % The model is for five degrees of freedom of the robot, 
  % with Foot1 touching the ground, 
  % q1: absolute angle of torso
  % q31: torso-femur 1; q32: torso - femur 2
  % q41: femur 1 - tibia 1; q42: femur 2 - tibia 2
  syms q31L q32L q41L q42L q1L real % relative angle
  syms dq31L dq32L dq41L dq42L dq1L real
  syms ddq31L ddq32L ddq41L ddq42L ddq1L
  syms g L1 L3 L4 M1 M3 M4 real
  syms MY1 MY3 MY4 MZ1 MZ3 MZ4 XX1 XX3 XX4 real
  syms rVPP gam
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
  % L3: Length of femur; L4: Length of Tibia
  R = @(t) [cos(t) -sin(t); sin(t) cos(t)]; % rotation 
  pH = -(R(q31)*[0;1]*L3 + R(q41)*[0;1]*L4);
  pG1 = pH + R(q31)*[0;1]*L3; % knee joint 1 
  pG2 = pH + R(q32)*[0;1]*L3; % knee joint 2
  %
  %
  vH =  jacobian(pH,q)*dq ;% 
  ddq = [ddq31L ddq32L ddq41L ddq42L ddq1L].';
  aH =  jacobian(vH,[q;dq])*[dq;ddq];  
  vG1 = jacobian(pG1,q)*dq;
  vG2 = jacobian(pG2,q)*dq;
  %
    % centers of gravity and positions of the various members 
  pT    = pH + R(q1)*[MY1; MZ1];
  pFem1 = pH + R(q31)*[MY3; MZ3];
  pFem2 = pH + R(q32)*[MY3; MZ3];
  pTib1 = pG1 + R(q41)*[MY4; MZ4];
  pTib2 = pG2 + R(q42)*[MY4; MZ4];
  pFoot1= pG1 + R(q41)*[0;1]*L4;
  pFoot2= pG2 + R(q42)*[0;1]*L4;
  pHead = pH +  R(q1)*[0;1]*L1;
  % velocity
  vFem1 = jacobian(pFem1,q)*dq;
  vFem2 = jacobian(pFem2,q)*dq;
  vTib1 = jacobian(pTib1,q)*dq;
  vTib2 = jacobian(pTib2,q)*dq;
  vFoot1=jacobian(pFoot1,q)*dq;
  vFoot2=jacobian(pFoot2,q)*dq;
  vT=jacobian(pT,q)*dq;
  %   
  KET    = simplify(1/2*M1*vT.'*vT ...
		     +1/2*XX1*dq1^2);
  KEFem1 = simplify(1/2*M3*vFem1.'*vFem1 ...
		    +1/2*XX3*dq31^2);
  KEFem2 = simplify(1/2*M3*vFem2.'*vFem2 ...
		    +1/2*XX3*dq32^2);		    
  KETib1 = simplify(1/2*M4*vTib1.'*vTib1 ...
		    +1/2*XX4*dq41^2);
  KETib2 = simplify(1/2*M4*vTib2.'*vTib2 ...
		    +1/2*XX4*dq42^2);
  %
  KE = simplify(KET + KEFem1 + KEFem2 + KETib1 + KETib2);
  % Potiental Energy
  PE = g*(M1*pT(2) + M3*pFem1(2) + M3*pFem2(2) + M4*pTib1(2) + ...
	  M4*pTib2(2));
  PE = simplify(PE);
  %
  %
  % Center of mass calculation...may want to use this for the theta1
  % parameter to watch the center of mass trajectory (see notes
  % 3/5/01 or Principles of Dynamics by Greenwood, page 136)
  %
  %
  M_total = M1 + 2*M3 + 2*M4;
  cm = 1/M_total*(M1*pT + M3*pFem1 + M3*pFem2 + M4*pTib1 + M4*pTib2);
  pCoM = cm; 
  pVPP = pCoM + R(gam)*[0;rVPP];
  cm = simplify(cm);
  dcm = jacobian(cm,q)*dq;
  acm = jacobian(dcm,[q;dq])*[dq;ddq];
  %
  % Calculate jacobian

  %
  %
  % Model NOTATION: Spong and Vidyasagar, page 142, Eq. (6.3.12)
  %                 D(q)ddq + C(q,dq)*dq + G(q) = B*tau + E2*F_external
  %
  %L=KE-PE;
  G = jacobian(PE,q).';
  G = simplify(G);
  D = jacobian(KE,dq).';
  D = simplify(D);
  D = jacobian(D,dq);
  D = simplify(D);
  
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
  %
  % Compute the matrix for the input torques and then the contact
  % points.
  %
  Phi_0 = [q31-q1;q32-q1;q41-q31;q42-q32;q1];
  B = jacobian(Phi_0,q);
  B = B.'*[eye(4,4);zeros(1,4)];
  %
  %
  % F_ext = [F_T^1;F_N^1;F_T^2;F_N^2];
  %
  Phi_1 = pFoot2;
  E1 = jacobian(Phi_1,q).';
  % Switching matrix --> at impact
  xi = [0  1  0  0  0;
        1  0  0  0  0;
        0  0  0  1  0;
        0  0  1  0  0;
        0  0  0  0  1];
  
  save pm_mat_files_sym/work_symb_model_5links_5DoF
  disp('[DONE saving MAT_FILES_SYM/work_symb_model_5links_5DoF.MAT]')
else
  load pm_mat_files_sym/work_symb_model_5links_5DOF
  disp('[DONE loading MAT_FILES_SYM/work_symb_model_5links_5DoF.MAT]')
end

N = max(size(q));

% First, output model to a file
%
% write file header
%
fid=print_function_header('dyn_mod_5links_5DoF',...
			  'D,C,G,B',...
			  'q,dq',...
			  ['is the kneed biped walking model.']);

%% Read in constants
%
fprintf(fid,'modelP;\n\n');

%% Reassign configuration parameters
%
fprintf(fid,'q31L=q(1); q32L=q(2); q41L=q(3); q42L=q(4); q1L=q(5);\n');
fprintf(fid,['dq31L=dq(1); dq32L=dq(2); dq41L=dq(3); dq42L=dq(4); dq1L=' ...
	     ' dq(5);\n\n']);

%% Model output
%
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


fclose(fid);


%% write file header
%
fid=print_function_header('swing_foot_height',...
			  'h',...
			  'x',...
			  ['']);

fprintf(fid,'modelP;\n\n');

fprintf(fid,'[nn,mm] = size(x);\n\n');
fprintf(fid,'if nn ~= 1 & mm ~= 1\n');
ttt = char(pFoot2(2));
ttt = rep_str(ttt,1);
fprintf(fid,'  h = %s;\n',ttt);
fprintf(fid,'else\n');
ttt = char(pFoot2(2));
ttt = rep_str(ttt,2);
fprintf(fid,'  h = %s;\n',ttt);
fprintf(fid,'end\n');
fclose(fid);


%% write file header
%
fid=print_function_header('head_height',...
			  'h',...
			  'x',...
			  ['']);

fprintf(fid,'modelP;\n\n');

fprintf(fid,'[nn,mm] = size(x);\n\n');
fprintf(fid,'if nn ~= 1 & mm ~= 1\n');
ttt = char(pHead(2));
ttt = rep_str(ttt,1);
fprintf(fid,'  h = %s;\n',ttt);
fprintf(fid,'else\n');
ttt = char(pHead(2));
ttt = rep_str(ttt,2);
fprintf(fid,'  h = %s;\n',ttt);
fprintf(fid,'end\n');
fclose(fid);


%% write file header
%
fid=print_function_header('limb_position',...
			  'out',...
			  'x,pH_horiz',...
			  ['']);

fprintf(fid,'modelP;\n\n');

pH_tmp = [0; pH(2)];

fprintf(fid,'[nn,mm] = size(x);\n');
fprintf(fid,'if nn == 10, k = mm; else, k = nn; end\n\n');
fprintf(fid,'if nn ~= 1 & mm ~= 1\n');
fprintf(fid,'  out.pH1 = pH_horiz*ones(k,1);\n');
ttt = char(pH(2));
ttt = rep_str(ttt,1);
fprintf(fid,'  out.pH2 = %s;\n',ttt);
pT_tmp = pT - pH;
ttt = char(pT_tmp(1));
ttt = rep_str(ttt,1);
fprintf(fid,'  out.pT1 = out.pH1 + %s;\n',ttt);
ttt = char(pT_tmp(2));
ttt = rep_str(ttt,1);
fprintf(fid,'  out.pT2 = out.pH2 + %s;\n',ttt);
% CoM
pCoM_tmp = pCoM - pH;
ttt = char(pCoM_tmp(1));
ttt = rep_str(ttt,1);
fprintf(fid,'  out.pCoM1 = out.pH1 + %s;\n',ttt);
ttt = char(pCoM_tmp(2));
ttt = rep_str(ttt,1);
fprintf(fid,'  out.pCoM2 = out.pH2 + %s;\n',ttt);
% VPP 
pVPP_tmp = pVPP - pH;
ttt = char(pVPP_tmp(1));
ttt = rep_str(ttt,1);
fprintf(fid,'  out.pVPP1 = out.pH1 + %s;\n',ttt);
ttt = char(pVPP_tmp(2));
ttt = rep_str(ttt,1);
fprintf(fid,'  out.pVPP2 = out.pH2 + %s;\n',ttt);
%
pHead_tmp = pHead - pH;
ttt = char(pHead_tmp(1));
ttt = rep_str(ttt,1);
fprintf(fid,'  out.pHead1 = out.pH1 + %s;\n',ttt);
ttt = char(pHead_tmp(2));
ttt = rep_str(ttt,1);
fprintf(fid,'  out.pHead2 = out.pH2 + %s;\n',ttt);
pFem1_tmp = pFem1 - pH;
ttt = char(pFem1_tmp(1));
ttt = rep_str(ttt,1);
fprintf(fid,'  out.pFem11 = out.pH1 + %s;\n',ttt);
ttt = char(pFem1_tmp(2));
ttt = rep_str(ttt,1);
fprintf(fid,'  out.pFem12 = out.pH2 + %s;\n',ttt);
pFem2_tmp = pFem2 - pH;
ttt = char(pFem2_tmp(1));
ttt = rep_str(ttt,1);
fprintf(fid,'  out.pFem21 = out.pH1 + %s;\n',ttt);
ttt = char(pFem2_tmp(2));
ttt = rep_str(ttt,1);
fprintf(fid,'  out.pFem22 = out.pH2 + %s;\n',ttt);
pG1_tmp = pG1 - pH;
ttt = char(pG1_tmp(1));
ttt = rep_str(ttt,1);
fprintf(fid,'  out.pG11 = out.pH1 + %s;\n',ttt);
ttt = char(pG1_tmp(2));
ttt = rep_str(ttt,1);
fprintf(fid,'  out.pG12 = out.pH2 + %s;\n',ttt);
pG2_tmp = pG2 - pH;
ttt = char(pG2_tmp(1));
ttt = rep_str(ttt,1);
fprintf(fid,'  out.pG21 = out.pH1 + %s;\n',ttt);
ttt = char(pG2_tmp(2));
ttt = rep_str(ttt,1);
fprintf(fid,'  out.pG22 = out.pH2 + %s;\n',ttt);
pTib1_tmp = pTib1 - pH;
ttt = char(pTib1_tmp(1));
ttt = rep_str(ttt,1);
fprintf(fid,'  out.pTib11 = out.pH1 + %s;\n',ttt);
ttt = char(pTib1_tmp(2));
ttt = rep_str(ttt,1);
fprintf(fid,'  out.pTib12 = out.pH2 + %s;\n',ttt);
pTib2_tmp = pTib2 - pH;
ttt = char(pTib2_tmp(1));
ttt = rep_str(ttt,1);
fprintf(fid,'  out.pTib21 = out.pH1 + %s;\n',ttt);
ttt = char(pTib2_tmp(2));
ttt = rep_str(ttt,1);
fprintf(fid,'  out.pTib22 = out.pH2 + %s;\n',ttt);
pFoot1_tmp = pFoot1 - pH;
ttt = char(pFoot1_tmp(1));
ttt = rep_str(ttt,1);
fprintf(fid,'  out.pFoot11 = out.pH1 + %s;\n',ttt);
ttt = char(pFoot1_tmp(2));
ttt = rep_str(ttt,1);
fprintf(fid,'  out.pFoot12 = out.pH2 + %s;\n',ttt);
pFoot2_tmp = pFoot2 - pH;
ttt = char(pFoot2_tmp(1));
ttt = rep_str(ttt,1);
fprintf(fid,'  out.pFoot21 = out.pH1 + %s;\n',ttt);
ttt = char(pFoot2_tmp(2));
ttt = rep_str(ttt,1);
fprintf(fid,'  out.pFoot22 = out.pH2 + %s;\n',ttt);
fprintf(fid,'else\n');
fprintf(fid,'  out.pH1 = pH_horiz;\n');
ttt = char(pH(2));
ttt = rep_str(ttt,2);
fprintf(fid,'  out.pH2 = %s;\n',ttt);
pT_tmp = pT - pH;
ttt = char(pT_tmp(1));
ttt = rep_str(ttt,2);
fprintf(fid,'  out.pT1 = out.pH1 + %s;\n',ttt);
ttt = char(pT_tmp(2));
ttt = rep_str(ttt,2);
fprintf(fid,'  out.pT2 = out.pH2 + %s;\n',ttt);
% CoM
pCoM_tmp = pCoM - pH;
ttt = char(pCoM_tmp(1));
ttt = rep_str(ttt,1);
fprintf(fid,'  out.pCoM1 = out.pH1 + %s;\n',ttt);
ttt = char(pCoM_tmp(2));
ttt = rep_str(ttt,1);
fprintf(fid,'  out.pCoM2 = out.pH2 + %s;\n',ttt);
% VPP 
pVPP_tmp = pVPP - pH;
ttt = char(pVPP_tmp(1));
ttt = rep_str(ttt,1);
fprintf(fid,'  out.pVPP1 = out.pH1 + %s;\n',ttt);
ttt = char(pVPP_tmp(2));
ttt = rep_str(ttt,1);
fprintf(fid,'  out.pVPP2 = out.pH2 + %s;\n',ttt);
%
pHead_tmp = pHead - pH;
ttt = char(pHead_tmp(1));
ttt = rep_str(ttt,2);
fprintf(fid,'  out.pHead1 = out.pH1 + %s;\n',ttt);
ttt = char(pHead_tmp(2));
ttt = rep_str(ttt,2);
fprintf(fid,'  out.pHead2 = out.pH2 + %s;\n',ttt);
pFem1_tmp = pFem1 - pH;
ttt = char(pFem1_tmp(1));
ttt = rep_str(ttt,2);
fprintf(fid,'  out.pFem11 = out.pH1 + %s;\n',ttt);
ttt = char(pFem1_tmp(2));
ttt = rep_str(ttt,2);
fprintf(fid,'  out.pFem12 = out.pH2 + %s;\n',ttt);
pFem2_tmp = pFem2 - pH;
ttt = char(pFem2_tmp(1));
ttt = rep_str(ttt,2);
fprintf(fid,'  out.pFem21 = out.pH1 + %s;\n',ttt);
ttt = char(pFem2_tmp(2));
ttt = rep_str(ttt,2);
fprintf(fid,'  out.pFem22 = out.pH2 + %s;\n',ttt);
pG1_tmp = pG1 - pH;
ttt = char(pG1_tmp(1));
ttt = rep_str(ttt,2);
fprintf(fid,'  out.pG11 = out.pH1 + %s;\n',ttt);
ttt = char(pG1_tmp(2));
ttt = rep_str(ttt,2);
fprintf(fid,'  out.pG12 = out.pH2 + %s;\n',ttt);
pG2_tmp = pG2 - pH;
ttt = char(pG2_tmp(1));
ttt = rep_str(ttt,2);
fprintf(fid,'  out.pG21 = out.pH1 + %s;\n',ttt);
ttt = char(pG2_tmp(2));
ttt = rep_str(ttt,2);
fprintf(fid,'  out.pG22 = out.pH2 + %s;\n',ttt);
pTib1_tmp = pTib1 - pH;
ttt = char(pTib1_tmp(1));
ttt = rep_str(ttt,2);
fprintf(fid,'  out.pTib11 = out.pH1 + %s;\n',ttt);
ttt = char(pTib1_tmp(2));
ttt = rep_str(ttt,2);
fprintf(fid,'  out.pTib12 = out.pH2 + %s;\n',ttt);
pTib2_tmp = pTib2 - pH;
ttt = char(pTib2_tmp(1));
ttt = rep_str(ttt,2);
fprintf(fid,'  out.pTib21 = out.pH1 + %s;\n',ttt);
ttt = char(pTib2_tmp(2));
ttt = rep_str(ttt,2);
fprintf(fid,'  out.pTib22 = out.pH2 + %s;\n',ttt);
pFoot1_tmp = pFoot1 - pH;
ttt = char(pFoot1_tmp(1));
ttt = rep_str(ttt,2);
fprintf(fid,'  out.pFoot11 = out.pH1 + %s;\n',ttt);
ttt = char(pFoot1_tmp(2));
ttt = rep_str(ttt,2);
fprintf(fid,'  out.pFoot12 = out.pH2 + %s;\n',ttt);
pFoot2_tmp = pFoot2 - pH;
ttt = char(pFoot2_tmp(1));
ttt = rep_str(ttt,2);
fprintf(fid,'  out.pFoot21 = out.pH1 + %s;\n',ttt);
ttt = char(pFoot2_tmp(2));
ttt = rep_str(ttt,2);
fprintf(fid,'  out.pFoot22 = out.pH2 + %s;\n',ttt);
fprintf(fid,'end\n');
fclose(fid);


%% write file header
%
fid=print_function_header('swing_foot_velocity',...
			  'vh,vv',...
			  'x',...
			  ['']);

fprintf(fid,'modelP;\n\n');

fprintf(fid,'[nn,mm] = size(x);\n\n');
fprintf(fid,'if nn ~= 1 & mm ~= 1\n');
ttt = char(vFoot2(1));
ttt = rep_str(ttt,1);
fprintf(fid,'  vh = %s;\n',ttt);
ttt = char(vFoot2(2));
ttt = rep_str(ttt,1);
fprintf(fid,'  vv = %s;\n',ttt);
fprintf(fid,'else\n');
ttt = char(vFoot2(1));
ttt = rep_str(ttt,2);
fprintf(fid,'  vh = %s;\n',ttt);
ttt = char(vFoot2(2));
ttt = rep_str(ttt,2);
fprintf(fid,'  vv = %s;\n',ttt);
fprintf(fid,'end');
fclose(fid);
%
%
%% write file header
%
fid=print_function_header('step_length',...
			  'h',...
			  'x',...
			  ['']);

fprintf(fid,'modelP;\n\n');
%
%% Reassign configuration parameters
%
fprintf(fid,'q31L=x(1); q32L=x(2); q41L=x(3); q42L=x(4); q1L=x(5);\n');

fprintf(fid,'h = %s;',fixlength(char(pFoot2(1)-pFoot1(1)),'+',70));
fclose(fid);
%
%
%% write file header
%
fid=print_function_header('hip_pos',...
			  'pHh,pHv',...
			  'x',...
			  ['']);

fprintf(fid,'modelP;\n\n');

fprintf(fid,'[nn,mm] = size(x);\n\n');
fprintf(fid,'if nn ~= 1 & mm ~= 1\n');
ttt = char(pH(1));
ttt = rep_str(ttt,1);
fprintf(fid,'  pHh = %s;\n',ttt);
ttt = char(pH(2));
ttt = rep_str(ttt,1);
fprintf(fid,'  pHv = %s;\n',ttt);
fprintf(fid,'else\n');
ttt = char(pH(1));
ttt = rep_str(ttt,2);
fprintf(fid,'  pHh = %s;\n',ttt);
ttt = char(pH(2));
ttt = rep_str(ttt,2);
fprintf(fid,'  pHv = %s;\n',ttt);
fprintf(fid,'end\n');
fclose(fid);
%
%
fid=print_function_header('CoM_pos',...
			  'pCh,pCv',...
			  'x',...
			  ['']);

fprintf(fid,'modelP;\n\n');

fprintf(fid,'[nn,mm] = size(x);\n\n');
fprintf(fid,'if nn ~= 1 & mm ~= 1\n');
ttt = char(pCoM(1));
ttt = rep_str(ttt,1);
fprintf(fid,'  pCh = %s;\n',ttt);
ttt = char(pCoM(2));
ttt = rep_str(ttt,1);
fprintf(fid,'  pCv = %s;\n',ttt);
fprintf(fid,'else\n');
ttt = char(pCoM(1));
ttt = rep_str(ttt,2);
fprintf(fid,'  pCh = %s;\n',ttt);
ttt = char(pCoM(2));
ttt = rep_str(ttt,2);
fprintf(fid,'  pCv = %s;\n',ttt);
fprintf(fid,'end\n');
fclose(fid);
%
%% write file header
%
fid=print_function_header('hip_vel',...
			  'vHh,vHv',...
			  'x',...
			  ['']);

fprintf(fid,'modelP;\n\n');

fprintf(fid,'[nn,mm] = size(x);\n\n');
fprintf(fid,'if nn ~= 1 & mm ~= 1\n');
ttt = char(vH(1));
ttt = rep_str(ttt,1);
fprintf(fid,'  vHh = %s;\n',ttt);
ttt = char(vH(2));
ttt = rep_str(ttt,1);
fprintf(fid,'  vHv = %s;\n',ttt);
fprintf(fid,'else\n');
ttt = char(vH(1));
ttt = rep_str(ttt,2);
fprintf(fid,'  vHh = %s;\n',ttt);
ttt = char(vH(2));
ttt = rep_str(ttt,2);
fprintf(fid,'  vHv = %s;\n',ttt);
fprintf(fid,'end\n');
fclose(fid);
%
%
%
%% write file header
%
fid=print_function_header('potential_energy_5DoF',...
			  'PE',...
			  'x',...
			  ['']);

fprintf(fid,'modelP;\n\n');

fprintf(fid,'[nn,mm] = size(x);\n\n');
fprintf(fid,'if nn ~= 1 & mm ~= 1\n');
ttt = char(PE);
ttt = rep_str(ttt,1);
fprintf(fid,'  PE = %s;\n',fixlength(ttt,'*+-',65,'       '));
fprintf(fid,'else\n');
ttt = char(PE);
ttt = rep_str(ttt,2);
fprintf(fid,'  PE = %s;\n',fixlength(ttt,'*+-',65,'       '));
fprintf(fid,'end\n');
fclose(fid);
%
%
%
%% write file header
%
fid=print_function_header('kinetic_energy_5DoF',...
			  'KE,PE',...
			  'x',...
			  ['']);

fprintf(fid,'modelP;\n\n');

fprintf(fid,'[nn,mm] = size(x);\n\n');
fprintf(fid,'if nn ~= 1 & mm ~= 1\n');
ttt = char(KE);
ttt = rep_str(ttt,1);
fprintf(fid,'  KE = %s;\n',fixlength(ttt,'*+-',65,'       '));
fprintf(fid,'else\n');
ttt = char(KE);
ttt = rep_str(ttt,2);
fprintf(fid,'  KE = %s;\n',fixlength(ttt,'*+-',65,'       '));
fprintf(fid,'end\n');
%
fprintf(fid,'[nn,mm] = size(x);\n\n');
fprintf(fid,'if nn ~= 1 & mm ~= 1\n');
ttt = char(PE);
ttt = rep_str(ttt,1);
fprintf(fid,'  PE = %s;\n',fixlength(ttt,'*+-',65,'       '));
fprintf(fid,'else\n');
ttt = char(PE);
ttt = rep_str(ttt,2);
fprintf(fid,'  PE = %s;\n',fixlength(ttt,'*+-',65,'       '));
fprintf(fid,'end\n');

%

fclose(fid);
%
%
%
%% write file header
%
fid=print_function_header('center_of_mass_5DoF',...
			  'cmh,cmv,dcmh,dcmv,acmh,acmv',...
			  'x,ddq',...
			  ['']);

fprintf(fid,'modelP;\n\n');

fprintf(fid,'[nn,mm] = size(x);\n\n');
fprintf(fid,'if nn ~= 1 & mm ~= 1\n');
ttt = char(cm(1));
ttt = rep_str(ttt,1);
fprintf(fid,'  cmh = %s;\n',fixlength(ttt,'*+-',65,'        '));
ttt = char(cm(2));
ttt = rep_str(ttt,1);
fprintf(fid,'  cmv = %s;\n',fixlength(ttt,'*+-',65,'        '));
ttt = char(dcm(1));
ttt = rep_str(ttt,1);
fprintf(fid,'  dcmh = %s;\n',fixlength(ttt,'*+-',65,'         '));
ttt = char(dcm(2));
ttt = rep_str(ttt,1);
fprintf(fid,'  dcmv = %s;\n',fixlength(ttt,'*+-',65,'         '));
ttt = char(acm(1));
ttt = rep_str(ttt,1);
fprintf(fid,'  if nargin > 1\n');
fprintf(fid,'    acmh = %s;\n',fixlength(ttt,'*+-',65,'           '));
ttt = char(acm(2));
ttt = rep_str(ttt,1);
fprintf(fid,'    acmv = %s;\n',fixlength(ttt,'*+-',65,'           '));
fprintf(fid,'  end\n');
fprintf(fid,'else\n');
ttt = char(cm(1));
ttt = rep_str(ttt,2);
fprintf(fid,'  cmh = %s;\n',fixlength(ttt,'*+-',65,'        '));
ttt = char(cm(2));
ttt = rep_str(ttt,2);
fprintf(fid,'  cmv = %s;\n',fixlength(ttt,'*+-',65,'        '));
ttt = char(dcm(1));
ttt = rep_str(ttt,2);
fprintf(fid,'  dcmh = %s;\n',fixlength(ttt,'*+-',65,'         '));
ttt = char(dcm(2));
ttt = rep_str(ttt,2);
fprintf(fid,'  dcmv = %s;\n',fixlength(ttt,'*+-',65,'         '));
ttt = char(acm(1));
ttt = rep_str(ttt,2);
fprintf(fid,'  if nargin > 1\n');
fprintf(fid,'    acmh = %s;\n',fixlength(ttt,'*+-',65,'           '));
ttt = char(acm(2));
ttt = rep_str(ttt,2);
fprintf(fid,'    acmv = %s;\n',fixlength(ttt,'*+-',65,'           '));
fprintf(fid,'  end\n');
fprintf(fid,'end\n');
fclose(fid);
%
%
%
fid=print_function_header('xi',...
			  'out',...
			  'in',...
			  ['']);
n=max(size(q));
fprintf(fid,'xi = zeros(%s);\n',num2str(n));
for k=1:n,
  for j=1:n,
    if xi(k,j)~=0
      ttt=num2str(xi(k,j));
      fprintf(fid,'xi(%s,%s) = %s;\n',num2str(k),num2str(j),ttt);
    end
  end
end
fprintf(fid,'out = xi*in;');
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
%%
