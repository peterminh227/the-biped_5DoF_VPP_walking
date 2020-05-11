function COT = walker_zd(alpha)
global optim_var
params = modelParameters;
L3 = params.L3; L4 = params.L4;
H0 = [1 0 0 0 0;0 1 0 0 0;0 0 1 0 0;0 0 0 1 0];
c = [-1,0,-1/2,0,-1];
H = [H0;c];
% calculate the dependent alpha
aa = [zeros(2,1);alpha(1:5)]; bb = [zeros(2,1);alpha(6:10)]; cc = [zeros(2,1);alpha(11:15)]; dd = [zeros(2,1);alpha(16:20)];
a_mat = [aa';bb';cc';dd'];
%% find qMinus
        function val = biped_configuration(x)
        temp = a_mat(:,end);
        q1 = temp(1); q2 = temp(2); q3 = temp(3); q4 = temp(4);
        % This is a geometric constraint to swing foot height = 0
        val =   L3*cos(q2 + x) - L3*cos(q1 + x) - L4*cos(q1 + q3 + x) + L4*cos(q2 + q4 + x);
        end
options = optimset('Display','Off');
x = fsolve(@biped_configuration,0,options);
qMinus = [a_mat(:,end);x];
thetaMinus = c*qMinus;
%% Get transition matrix
S = [0  1  0  0  0;
     1  0  0  0  0;
     0  0  0  1  0;
     0  0  1  0  0;
     0  0  0  0  1];
qPlus = S*qMinus;
[De,~,~,~,E] = dyn_mod_5links_7DoF([qMinus;zeros(2,1)],zeros(1:7));
E = E(:,3:4)';
Lambda = eye(7) - inv(De)*E'*inv(E*inv(De)*E')*E;
[~,~,~,~,E1] = dyn_mod_5links_7DoF([qPlus;zeros(2,1)],zeros(1:7));
E1 = E1(:,3:4)';  % take the first row
A = S*(Lambda(1:5,1:5)-Lambda(1:5,6:7)*E1(:,1:5)*S); % transition matrix from dx- --> dx+
%% calculate thetaPlus
temp = H*S*inv(H)*[a_mat(:,end); thetaMinus];
alpha_x(:,1) = temp(1:4);
thetaPlus = temp(5);
m = thetaMinus - thetaPlus;
% calculate alpha_1
M = size(a_mat,2)-1;
omega = H\[M*(a_mat(:,end)-a_mat(:,end-1))/(thetaMinus-thetaPlus); 1];
alpha_x(:,2) = (thetaMinus-thetaPlus)*H0*A*omega/(M*c*A*omega) + alpha_x(:,1);
a_mat(:,1:2) = alpha_x(:,1:2);
alpha = [a_mat(1,:)'; a_mat(2,:)';a_mat(3,:)';a_mat(4,:)']; % full of Bezier coefficient
%% -> till this step, we have the dependent coefficient of two first column 
% calculate fixed point xi2 -> important value
D_plus = dyn_mod_5links_5DoF(qPlus,zeros(5,1)); 
D_minus = dyn_mod_5links_5DoF(qMinus,zeros(5,1));
D_Nplus = D_plus(5,:); D_Nminus = D_minus(5,:);
jac_h = jacobianH(qMinus,zeros(5,1),alpha,thetaPlus,m);
velTrans = [jac_h;D_Nminus];
lambdaBar = velTrans\[zeros(5-1,1);1];
deltaZero = D_Nplus*A*lambdaBar;
C1 = deltaZero ;%> stability value
xi1 = linspace(thetaPlus,thetaMinus,1e5); % thetaPlus->thetaPlus
alpha_times_beta = alpha_beta(xi1,alpha,thetaPlus,m); %kappa1/kappa2
C2_cum = cumtrapz(xi1,alpha_times_beta);
K = max(C2_cum); % use for non-linear condition
C2 = trapz(xi1,alpha_times_beta);
xi2_minus = -sqrt(2*C2/(1-C1^2)); % fixed point
optim_var = [alpha;thetaPlus;m;xi2_minus];
% condition 1/2*xi2_plus^2*C1^2-K >0
%% --> we have xi2_minus, x12_plus
[~,~,~,~,~,~,COT,~] = sim_zd;
if (COT <0)
    COT = abs(COT)*1e10;
end
end