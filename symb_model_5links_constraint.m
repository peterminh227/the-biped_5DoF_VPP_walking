function symb_model_5links_constraint
    % run symb_model_5links_7DoF first
    load pm_mat_files_sym/work_symb_model_5links_7DoF
    c1 = pFoot1;
    c2 = pFoot2;
    jac_c1 = jacobian(c1,q);         % c_dot  = jac_c * dq
    jac_c2 = jacobian(c2,q);
    gamma_c1 = jacobian(jac_c1*dq,q)*dq; % c_ddot = gamma_c + jac_c * ddq
    gamma_c2 = jacobian(jac_c2*dq,q)*dq;
    
    %
fid=print_function_header('constraint_5links_7DoF',...
			  'jac_c1,jac_c2,gamma_c1,gamma_c2,pFoot1,pFoot2,KE,PE',...
			  'q,dq',...
			  ['return the constrain value.']);
fprintf(fid,'modelP;\n\n');
fprintf(fid,['\nq31L=q(1);q32L=q(2);q41L=q(3);q42L=q(4);q1L=q(5);' ...
	       'ycm=q(6);zcm=q(7);']);
fprintf(fid,['\ndq31L=dq(1);dq32L=dq(2);dq41L=dq(3);dq42L=dq(4);' ...
	       'dq1L=dq(5);dycm=dq(6);dzcm=dq(7);']);
%%%%%
% jac_1
[n,m]=size(jac_c1);
fprintf(fid,'\n%s',['jac_c1=zeros(',num2str(n),',',num2str(m),');']);
for i=1:n
  for j=1:m
    Temp0=jac_c1(i,j);
    if Temp0~=0
      Temp1=char(Temp0);
      Temp2=['jac_c1(',num2str(i),',',num2str(j),')=',Temp1,';'];
      fprintf(fid,'\n%s',Temp2);
    end
  end
end
% jac_2
[n,m]=size(jac_c2);
fprintf(fid,'\n%s',['jac_c2=zeros(',num2str(n),',',num2str(m),');']);
for i=1:n
  for j=1:m
    Temp0=jac_c2(i,j);
    if Temp0~=0
      Temp1=char(Temp0);
      Temp2=['jac_c2(',num2str(i),',',num2str(j),')=',Temp1,';'];
      fprintf(fid,'\n%s',Temp2);
    end
  end
end
%
[n,m]=size(gamma_c1);
fprintf(fid,'\n%s',['gamma_c1=zeros(',num2str(n),',',num2str(m),');']);
for i=1:n
  for j=1:m
    Temp0=gamma_c1(i,j);
    if Temp0~=0
      Temp1=char(Temp0);
      Temp2=['gamma_c1(',num2str(i),',',num2str(j),')=',Temp1,';'];
      fprintf(fid,'\n%s',Temp2);
    end
  end
end
%
[n,m]=size(gamma_c2);
fprintf(fid,'\n%s',['gamma_c2=zeros(',num2str(n),',',num2str(m),');']);
for i=1:n
  for j=1:m
    Temp0=gamma_c2(i,j);
    if Temp0~=0
      Temp1=char(Temp0);
      Temp2=['gamma_c2(',num2str(i),',',num2str(j),')=',Temp1,';'];
      fprintf(fid,'\n%s',Temp2);
    end
  end
end
%
[n,m]=size(pFoot1);
fprintf(fid,'\n%s',['pFoot1=zeros(',num2str(n),',',num2str(m),');']);
for i=1:n
  for j=1:m
    Temp0=pFoot1(i,j);
    if Temp0~=0
      Temp1=char(Temp0);
      Temp2=['pFoot1(',num2str(i),',',num2str(j),')=',Temp1,';'];
      fprintf(fid,'\n%s',Temp2);
    end
  end
end
%
[n,m]=size(pFoot2);
fprintf(fid,'\n%s',['pFoot2=zeros(',num2str(n),',',num2str(m),');']);
for i=1:n
  for j=1:m
    Temp0=pFoot2(i,j);
    if Temp0~=0
      Temp1=char(Temp0);
      Temp2=['pFoot2(',num2str(i),',',num2str(j),')=',Temp1,';'];
      fprintf(fid,'\n%s',Temp2);
    end
  end
end
%
[n,m]=size(PE);
fprintf(fid,'\n%s',['PE=zeros(',num2str(n),',',num2str(m),');']);
for i=1:n
  for j=1:m
    Temp0=PE(i,j);
    if Temp0~=0
      Temp1=char(Temp0);
      Temp2=['PE(',num2str(i),',',num2str(j),')=',Temp1,';'];
      fprintf(fid,'\n%s',Temp2);
    end
  end
end
%
[n,m]=size(KE);
fprintf(fid,'\n%s',['KE=zeros(',num2str(n),',',num2str(m),');']);
for i=1:n
  for j=1:m
    Temp0=KE(i,j);
    if Temp0~=0
      Temp1=char(Temp0);
      Temp2=['KE(',num2str(i),',',num2str(j),')=',Temp1,';'];
      fprintf(fid,'\n%s',Temp2);
    end
  end
end
fclose(fid); 



end
