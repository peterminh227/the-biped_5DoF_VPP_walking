function symb_stance_force
%SYMB_STANCE_FORCE.M

%Copyright (c) 2003 by Eric Westervelt and Jessy Grizzle.  This code
%may be freely used for noncommercial ends. If use of this code in
%part or in whole results in publication, proper citation must be
%included in that publication.  This code comes with no guarantees or
%support.

%Eric R. Westervelt & Jessy W. Grizzle
%7/12/01
%11/30/03 - updated for the CDC

if isempty(dir('pm_mat_files_sym/work_symb_stance_force.mat'))
  % This is for a robot with an upright trunk, two legs and
  % knees. The model is for five degrees of freedom of the robot,
  % with Foot1 touching the ground
  %
  syms q31L q32L q41L q42L q1L real
  syms dq31L dq32L dq41L dq42L dq1L real
  syms ddq31L ddq32L ddq41L ddq42L ddq1L real
  syms y z dy dz ddy ddz real
  syms g L1 L3 L4 M1 M3 M4 real
  syms MY1 MY3 MY4 MZ1 MZ3 MZ4 XX1 XX3 XX4 IA3 IA4 real
  
  syms Lb1 Lt XXb XXc XXw PEbc PEcc PEwc real
  %
  % Change coordinates
  %
  [q_rel,dq_rel,q_abs,dq_abs] = rel2abs;
  %
  % old absolute coordinates in terms of lapin coordinates
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
  %
  % reference = absolute angles for everything; trigonometric sense for all
  %              angles, measured from the vertical,
  %
  q = [q; y; z];
  dq = [dq; dy; dz];
  ddq = [ddq31L; ddq32L; ddq41L; ddq42L; ddq1L; ddy; ddz];
  %
  % 
  % The rotation matrix, R(t), rotates the vector [MYx; MZx]/Mx (which
  % originates at the joint and ends at the link center of mass) an
  % angle t in the (Y-Z) plane (X points our of the plane).
  %
  R = @(t) [cos(t) -sin(t); sin(t) cos(t)]; % rotation 
  pH = [y;z]-(R(q31)*[0;1]*L3 + R(q41)*[0;1]*L4);
  pG1 = pH + R(q31)*[0;1]*L3;
  pG2 = pH + R(q32)*[0;1]*L3;
  %
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
  %
  % Model NOTATION: Spong and Vidyasagar, page 142, Eq. (6.3.12)
  %                 D(q)ddq + C(q,dq)*dq + G(q) = B*tau + E2*F_external
  %
  %L = KE-PE;
  G = jacobian(PE,q).';
  G = simplify(G);
  D = simplify(jacobian(KE,dq).');
  D = simplify(jacobian(D,dq));
  %
  syms C real
  n = max(size(q));
  for k = 1:n
    for j = 1:n
      C(k,j) = 0*g;
      for i = 1:n
	C(k,j) = C(k,j)+1/2*(diff(D(k,j),q(i)) + ...
			     diff(D(k,i),q(j))- ...
			     diff(D(i,j),q(k)))*dq(i);
      end
    end
  end
  C = simplify(C);
  %
  for k = 1:n
    for j = 1:n
      for i = 1:n
	cc(i,j,k) = 1/2*(diff(D(k,j),q(i)) + ...
			 diff(D(k,i),q(j)) - ...
			 diff(D(i,j),q(k)));
      end
    end
  end
  %
  % Compute the matrix for the input torques and then the contact
  % points.
  %
  Phi_0 = [q31-q1;q32-q1;q41-q31;q42-q32;q1;y;z];
  B = jacobian(Phi_0,q);
  B = B.'*[eye(4,4);zeros(3,4)];
  %
  % Compute matrices associated with the viscous and static
  % friction
  %
  D21 = D(6:7,1:5);
  C21 = C(6:7,1:5);
  G2  = G(6:7);
  %
  % See my notes, 7/13/01 for equations...
  %
  F = simplify(D21*ddq(1:5) + C21*dq(1:5) + G2);
  %
  save pm_mat_files_sym/work_symb_stance_force
else
  load pm_mat_files_sym/work_symb_stance_force
  disp('[DONE loading mat_files_sym/work_symb_stance_force.mat]')
end
%
%
%
%% write file header
%
fid=print_function_header('stance_force',...
			  'f_tan,f_norm',...
			  'x,ddq',...
			  ['']);

%% Read in constants
%
fprintf(fid,'modelP;\n\n');

fprintf(fid,'%% tangential force\n');
ttt = char(F(1));
ttt = rep_str(ttt,2);
ttt = fixlength(ttt,'*+-',65,'        ');
fprintf(fid,'f_tan = %s;\n\n',ttt);

fprintf(fid,'%% normal force\n');
ttt = char(F(2));
ttt = rep_str(ttt,2);
ttt = fixlength(ttt,'*+-',65,'         ');
fprintf(fid,'f_norm = %s;',ttt);
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
