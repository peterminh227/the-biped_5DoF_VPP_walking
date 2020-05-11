% OPTIMIZATION
clc;
global optim_var

% starting with guess value, v_step = 0.8 m/s
% Bezier degree 6; 
%aa = [  3.497673; 3.269202; 3.068850; 3.028060; 2.979016];
%bb = [  3.033769; 3.363808; 3.537780; 3.602939; 3.577896];
%cc = [ -0.334434; -0.223849; -0.138869; -0.119916; -0.143030];
%dd = [ -0.045280; -1.134701; -0.324863; -0.346557; -0.322206];
aa = [  3.5; 3.2; 3; 3; 2.98];
bb = [  3; 3.4; 3.5; 3.6; 3.58];
cc = [ -0.34; -0.22; -0.14; -0.1; -0.14];
dd = [ -0.05; -1.13; -0.32; -0.34; -0.32];
alpha_input0=[aa;bb;cc;dd];
%
nState = size(alpha_input0,1);
Aineq = []; Bineq = []; Aeq = []; Beq = []; 
LB = -2*pi*ones(20,1);
UB = 2*pi*ones(20,1);
LB(5) = 3*pi/4; UB(5) = 5*pi/4;
LB(10) = 3*pi/4; UB(10) = 5*pi/4;
LB(15) = -pi/2; UB(15) = 0;
LB(20) = -pi/2; UB(20) = 0;
NLP_display = 'iter-detailed';
Accuracy_type = 'low';
    switch Accuracy_type
        case 'low'
            options = optimset(...
                'Display',NLP_display,...
                'TolFun',1e-4,...
                'MaxIter',100,...
                'Algorithm','sqp',...
                'MaxFunEvals',1e4*(nState));
                
        case 'medium'
            options = optimset(...
                'Display',NLP_display,...
                'TolFun',1e-6,...
                'MaxIter',20000,...
                'MaxFunEvals',5e4*(nState));
        case 'high'
            options = optimset(...
                'Display',NLP_display,...
                'TolFun',1e-8,...
                'MaxIter',20000,...
                'MaxFunEvals',1e5*(nState));
        otherwise
            error('Invalid value for options.defaultAccuracy')
    end
%[presult,optfval] = fmincon(@walker_zd,alpha_input0,Aineq,Bineq,Aeq,Beq,LB,UB,@walker_cond,options);
%% print result
print_optimization(optim_var);
%% running parameter from optimization
fprintf('******* Collect data ******* \n');
%[a,b,m,gamma] = optimizationParameters;
%optim_var = [a;b;m;gamma];
[t,z,step_t,te,ze,ie,J,step_length] = sim_zd_af;
I = 0;
nI = length(te);
for k = 1:nI, I(k+1) = find(t == te(k)); end
y = [];
dy = [];
f = [];
x = [];
u = [];
for j = 1:nI
  x(I(j)+1:I(j+1),:) = zd_s_lift(z(I(j)+1:I(j+1),1:2),optim_var);
  for k = I(j)+1:I(j+1)
    [h,lfh] = decouple_5links_5DoF(x(k,1:5)',x(k,6:10)',optim_var);
    y = [y; h.'];
    dy = [dy; lfh.'];
    dx = walk(t(k),x(k,:)');
    [f_tan,f_norm] = stance_force(x(k,1:10),dx(1:10));
    f = [f;f_tan f_norm];
    u_tmp = control_synth(x(k,:),optim_var);
    u = [u; u_tmp];
  end
end
RUN_NAME = 'optim_traj';
save(['pm_mat_files/',RUN_NAME,'_zd_data'], ...
     't','z','x','u','y','dy','f','step_t');

doplots = [1; % positions
           0; % velocities
	  0; % torques
	   0; % outputs
	   1; % foot and hip trajectories
	   0; % hip specs & knee deflection specs
	   1]; % stance let forces

doprints = 0;

make_plots(t,x,u,y,dy,f,step_t,doplots,doprints,'ps');

anim(t,x,0.01);
