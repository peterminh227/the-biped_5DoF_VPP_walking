

% 11/2016 Peter Minh -- rigig body for the VPP model -  5 links
clc
global optim_var param 
% for the VPP controller

param.alpha0 = deg2rad(110);
param.L0 = 0.36;
param.k = 5000*1.5;
param.c = 10;
param.offset = -pi/40;
param.numSteps = 5;
%
%[a,b,m,gamma] = optimizationParameters;
%optim_var = [a;b;m;gamma];
RUN_NAME = 'test';
disp(['RUN_NAME = ',RUN_NAME]);
res = [];
% select controller
HZD = 0;
y = [];
dy = [];
if (HZD  == 1 ) 
    [t,x,u,y,dy,f,step_t] = walk_HZD;
else
    [t,x,u,f,step_t,VBLAout] = walk_VPP;
end
save(['pm_mat_files/',RUN_NAME,'_data'], ...
     't','x','u','y','dy','f','step_t');
doplots = [0; % positions
           1; % velocities
           1; % torques
           0; % outputs
           0; % foot and hip trajectories
           1; % hip specs & knee deflection specs
           1]; % stance let forces
x_red = [x(:,1:5),x(:,8:12)];
doprints = 0;
make_plots(t,x_red,u,y,dy,f,step_t,doplots,doprints,'fig');
% calculate ground reaction force
anim(t,x_red,0.005);