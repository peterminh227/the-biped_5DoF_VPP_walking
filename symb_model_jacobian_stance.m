function symb_model_jacobian_stance
    
    load pm_mat_files_sym/work_symb_model_5links_7DoF
    ST_p_B_TS = pFoot1 - pH;
    ST_p_T_S = R(q31L+q1L-pi)\ST_p_B_TS;
    ST_theta = atan2(ST_p_T_S(2),ST_p_T_S(1)) + q31L -pi+pi/2;
    ST_length_T_S = simplify(norm(ST_p_T_S));
    ST_J_polar_T_S = simplify(jacobian([ST_length_T_S;ST_theta],q),5);
    % for swing leg
    SW_p_B_TS = pFoot2 - pH;
    SW_p_T_S = R(q32L+q1L-pi)\SW_p_B_TS;
    SW_theta = atan2(SW_p_T_S(2),SW_p_T_S(1)) + q32L -pi+pi/2;
    SW_length_T_S = simplify(norm(SW_p_T_S));
    SW_J_polar_T_S = simplify(jacobian([SW_length_T_S;SW_theta],q),5);
%
fid=print_function_header('jacobian_stance',...
			  'ST_length_T_S,ST_theta,ST_J_polar_T_S,SW_length_T_S,SW_theta,SW_J_polar_T_S',...
			  'q,dq',...
			  ['mapping the BTSLIP to the rigid body']);

%% Read in constants
%
fprintf(fid,'modelP;\n\n');

%% Reassign configuration parameters
%
fprintf(fid,'q31L=q(1); q32L=q(2); q41L=q(3); q42L=q(4); q1L=q(5);\n');
fprintf(fid,['dq31L=dq(1); dq32L=dq(2); dq41L=dq(3); dq42L=dq(4); dq1L=' ...
	     ' dq(5);\n\n']);

%% Model output
fprintf(fid,'\n%% ST_theta, ST_length_T_S\n');
ttt=char(ST_length_T_S);
fprintf(fid,'ST_length_T_S = %s;\n',ttt);
ttt=char(ST_theta);
fprintf(fid,'ST_theta = %s;\n',ttt);
%
fprintf(fid,'\n%% ST_J_polar_T_S\n');
[N,M]=size(ST_J_polar_T_S);
fprintf(fid,'ST_J_polar_T_S=zeros(%s,%s);\n',num2str(N),num2str(M));
for k=1:N,
  for j=1:M,
    if ST_J_polar_T_S(k,j)~=0
      ttt=char(ST_J_polar_T_S(k,j));
      fprintf(fid,'ST_J_polar_T_S(%s,%s)=%s;\n',num2str(k),num2str(j),ttt);
    end
  end
end
%
fprintf(fid,'\n%% SW_theta, SW_length_T_S\n');
ttt=char(SW_length_T_S);
fprintf(fid,'SW_length_T_S = %s;\n',ttt);
ttt=char(SW_theta);
fprintf(fid,'SW_theta = %s;\n',ttt);
%
fprintf(fid,'\n%% SW_J_polar_T_S\n');
[N,M]=size(SW_J_polar_T_S);
fprintf(fid,'SW_J_polar_T_S=zeros(%s,%s);\n',num2str(N),num2str(M));
for k=1:N,
  for j=1:M,
    if SW_J_polar_T_S(k,j)~=0
      ttt=char(SW_J_polar_T_S(k,j));
      fprintf(fid,'SW_J_polar_T_S(%s,%s)=%s;\n',num2str(k),num2str(j),ttt);
    end
  end
end

fclose(fid);

end
