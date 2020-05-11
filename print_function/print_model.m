%PRINT_MODEL.M

%Copyright (c) 2003 by Eric Westervelt and Jessy Grizzle.  This code
%may be freely used for noncommercial ends. If use of this code in
%part or in whole results in publication, proper citation must be
%included in that publication.  This code comes with no guarantees or
%support.

%Eric R. Westervelt & Jessy W. Grizzle
%11/30/03 - updated for the CDC

fid=print_function_header(fcn_name,...
			  'D,C,G,B,E1,EFV,EFS',...
			  'q,dq',...
			  ['Model NOTATION: Spong and Vidyasagar, page ',...
		    '142, Eq. (6.3.12)     D(q)ddq + C(q,dq)*dq + G(q) =',...
		    ' B*tau + E*F_external']);

fprintf(fid,'modelP;\n\n');
if n == 5
  fprintf(fid,'q31L=q(1);q32L=q(2);q41L=q(3);q42L=q(4);q1L=q(5);');
  fprintf(fid,['\ndq31L=dq(1);dq32L=dq(2);dq41L=dq(3);dq42L=' ...
	       'dq(4);dq1L=dq(5);']);
else
  fprintf(fid,['\nq31L=q(1);q32L=q(2);q41L=q(3);q42L=q(4);q1L=q(5);' ...
	       'ycm=q(6);zcm=q(7);']);
  fprintf(fid,['\ndq31L=dq(1);dq32L=dq(2);dq41L=dq(3);dq42L=dq(4);' ...
	       'dq1L=dq(5);dycm=dq(6);dzcm=dq(7);']);
end
N = max(size(q));
fprintf(fid,'%% D matrix\n');
fprintf(fid,'D=zeros(%s);\n',num2str(N));
for i=1:N
  for j=1:N
    %{
    Temp0=D(i,j);
    if Temp0 ~= 0
      [r,s] = subexpr_mood(Temp0);
      fprintf(fid,'\nclear s\n');
      for k = 1:length(s)
	fprintf(fid,'s(%s) = %s;\n',num2str(k), ...
		fixlength(char(s(k)),'*+-',63,'       '));
      end
      Temp1=fixlength(char(r),'*+-',63,'       ');
      Temp2=['D(',num2str(i),',',num2str(j),')=',Temp1,';'];
      fprintf(fid,'\n%s',Temp2);
    end
    %}
    if D(i,j)~=0
      ttt=fixlength(char(D(i,j)),'*+-',65,'       ');
      fprintf(fid,'D(%s,%s)=%s;\n',num2str(i),num2str(j),ttt);
    end
  end
end
fprintf(fid,'\n%% C matrix\n');
fprintf(fid,'C=zeros(%s);\n',num2str(N));
for i=1:n
  for j=1:n
    %{
    Temp0=C(i,j);
    if Temp0 ~= 0
      [r,s] = subexpr_mood(Temp0);
      fprintf(fid,'\nclear s\n');
      for k = 1:length(s)
	fprintf(fid,'s(%s) = %s;\n',num2str(k), ...
		fixlength(char(s(k)),'*+-',63,'       '));
      end
      Temp1=fixlength(char(r),'*+-',63,'       ');
      Temp2=['C(',num2str(i),',',num2str(j),')=',Temp1,';'];
      fprintf(fid,'\n%s',Temp2);
    end
      %}
    if C(i,j)~=0
      ttt=fixlength(char(C(i,j)),'*+-',65,'       ');    
      fprintf(fid,'C(%s,%s)=%s;\n',num2str(i),num2str(j),ttt);
    end
  end
end
fprintf(fid,'\n%% G matrix\n');
fprintf(fid,'G=zeros(%s,1);\n',num2str(N));
for i=1:n
  %{
  [r,s] = subexpr_mood(G(i));
  if r ~= 0
    if length(s) > 0
      fprintf(fid,'\nclear s\n');
      for k = 1:length(s)
	fprintf(fid,'s(%s) = %s;\n',num2str(k), ...
		fixlength(char(s(k)),'*+-',63,'       '));
      end
    end
    Temp1=fixlength(char(r),'*+-',63,'       ');
  else
    Temp1 = '0';
  end
  Temp2=['G(',num2str(i),')=',Temp1,';'];
  fprintf(fid,'\n%s',Temp2);
    %}
  if G(i)~=0
    ttt=fixlength(char(G(i)),'*+-',65,'     ');
    fprintf(fid,'G(%s)=%s;\n',num2str(i),ttt);
  end
end
fprintf(fid,'\n%%');
fprintf(fid,'\n%%');
fprintf(fid,'\n%%');
[n,m]=size(B);
fprintf(fid,'\n%s',['B=zeros(',num2str(n),',',num2str(m),');']);
for i=1:n
  for j=1:m
    Temp0=B(i,j);
    if Temp0 ~= 0
      Temp1=char(Temp0);
      Temp2=['B(',num2str(i),',',num2str(j),')=',Temp1,';'];
      fprintf(fid,'\n%s',Temp2);
    end
  end
end
fprintf(fid,'\n%%');
fprintf(fid,'\n%%');
fprintf(fid,'\n%%');
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
fprintf(fid,'\n%%');
fprintf(fid,'\n%%');
fprintf(fid,'\n%%');
[n,m]=size(EFV);
fprintf(fid,'\n%s',['EFV=zeros(',num2str(n),',',num2str(m),');']);
for i=1:n
  Temp0=EFV(i);
  if Temp0~=0
    Temp1=char(Temp0);
    Temp2=['EFV(',num2str(i),')=',Temp1,';'];
    fprintf(fid,'\n%s',Temp2);
  end
end
fprintf(fid,'\n%%');
fprintf(fid,'\n%%');
fprintf(fid,'\n%%');
[n,m]=size(EFS);
fprintf(fid,'\n%s',['EFS=zeros(',num2str(n),',',num2str(m),');']);
for i=1:n
  Temp0=EFS(i);
  if Temp0~=0
    Temp1=char(Temp0);
    Temp2=['EFS(',num2str(i),')=',Temp1,';'];
    fprintf(fid,'\n%s',Temp2);
  end
end
fprintf(fid,'\n');
status = fclose(fid);
