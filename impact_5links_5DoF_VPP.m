function out = impact_5links_5DoF_VPP(x)
%IMPACT_LAP_RED    Map pre-impact state to post-impact state.
%

% 5 DoF
q = x(1:5);
dq = x(8:12);
% 
[num_r,num_c] = size(q);
if num_r < num_c
  q = q';
  dq = dq';
end
if (size(x,1) == 1)
    x_f = x';
else x_f = x;
end

qe = x_f(1:7);
dqe_minus = x_f(8:14);

[De,~,~,~,E1e] = dyn_mod_5links_7DoF(qe,dqe_minus);

E = E1e(:,3:4);
% De*(dq^+ - dq^-) = E*impulsive forces=E*Fe
% E'*dq^+ = 0;
%
A = [De,-E;E',zeros(2,2)];
B = [De*dqe_minus; zeros(2,1)];
temp = inv(A)*B;

% to test that after impact the fixed foot raises itself from the
% ground

[vh,vv] = stance_foot_vel([x_f(1:7); temp(1:7)]);

Fe = temp(8:9);
test1 = min(100,abs(Fe(1)/Fe(2)));
test2 = vv;
%if test1 <= 2/3 && test2 >= 0
  flag = 0;
  q_plus_new = xi(q);
  dq_plus_new = xi(temp(1:5));
  
  temp_a = red2full([q_plus_new;dq_plus_new]);
  %temp_a(6) = x(6);
  %temp_a(13) =x(13);
  %temp_a(14) =x(14);
  out = [temp_a;flag;test1;test2];

%else
%  error('swing leg end ground contact violation');
%end
