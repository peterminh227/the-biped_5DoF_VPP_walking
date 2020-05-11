function out = modelParameters
% parameter from JW
out.g = 9.81;

out.L1=0.4;
out.L3=0.19;
out.L4=0.19;

out.M1 = 12.5;
out.M3= 0.7;
out.M4= 0.7;
out.MY1=0;
out.MY3=0;
out.MY4=0;
out.MZ1= 0.21;
out.MZ3= 0.12577;
out.MZ4= 125.77/1000; 
out.XX1 = out.M1/12*(0.42^2+0.2^2); % inertial of torso
out.XX3= 447.42e-5; % femur
out.XX4= 447.42e-5;  % tibia
% VPP positon
out.rVPP = 0.1;
out.gam = pi/3;
% from the BTSLIP model
%{ 
out.g = 9.81;
out.L1=0.8;
out.L3=0.5;
out.L4=0.5;

out.M1=60;
out.M3=5;
out.M4=5;
out.MY1=0;
out.MY3=0;
out.MY4=0;
out.MZ1=0.35;
out.MZ3=0.25;
out.MZ4=0.25; 
out.XX1=1; % inertial of torso
out.XX3=0.1; % femur
out.XX4=0.1;  % tibia
% VPP positon
out.rVPP = 0.1;
out.gam = pi/3;
%}
%{

out.L1=0.63;
out.L3=0.4;
out.L4=0.4;

out.M1=12;
out.M3=6.8;
out.M4=3.2;

out.MY1=0;
out.MY3=0;
out.MY4=0;
out.MZ1=0.24;
out.MZ3=0.11;
out.MZ4=0.24; 
out.XX1=1.33; % inertial of torso
out.XX3=0.47; % femur
out.XX4=0.2;  % tibia
%}



