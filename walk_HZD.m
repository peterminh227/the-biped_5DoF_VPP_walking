function varargout = walk_HZD(t,x,flag)
%WALK   Simulate the kneed biped walker.

%Copyright (c) 2003 by Eric Westervelt and Jessy Grizzle.  This code
%may be freely used for noncommercial ends. If use of this code in
%part or in whole results in publication, proper citation must be
%included in that publication.  This code comes with no guarantees or
%support.

%Eric R. Westervelt & Jessy W. Grizzle
%10/2/00
%11/30/03 - updated for the CDC

if nargin == 2, flag = ''; end

if nargin == 0, flag = 'demo'; end

switch flag
 case ''                                 % Return dx/dt = dynamics(t,x).
  varargout{1} = f(t,x);
 case 'events'                           % Return [value,isterminal,direction].
  [varargout{1:3}] = events(t,x);
 case 'demo'                             % Run a demo.
  [varargout{1:7}] = demo;
 otherwise
  error(['Unknown flag ''' flag '''.']);
end


% --------------------------------------------------------------------------
% This is the system dynamics function

function dx = f(~,x)
global param
global optim_var
q = x(1:5); dq = x(8:12);
q_f = x(1:7);
dq_f = x(8:14);
HZD =1; % HZD controller

[H,LfH,L2fH,LgLfH] = decouple_5links_5DoF(q,dq,optim_var);
[D_f,C_f,G_f,B_f] = dyn_mod_5links_7DoF(q_f,dq_f);
[jac_c1,jac_c2,gamma_c1,gamma_c2,~,~] = constraint_5links_7DoF(q_f,dq_f);
%[KE,PE] = kinetic_energy_5DoF([q',dq']);
gains = [1, 1, 1, 1];
out = [bb([H(1),LfH(1)],gains(1));
       bb([H(2),LfH(2)],gains(2));
       bb([H(3),LfH(3)],gains(3));
       bb([H(4),LfH(4)],gains(4))];
psi = out(:,1);

% Used for controller that use feedback linearization
u = LgLfH\(psi-L2fH);


%ddq = inv_D*(-C*dq-G+B*u);
% D*ddq + C*dq + G = B*u +J'*F
FSM = 1;

if FSM ==1 % Single support
    Af = [D_f,-jac_c1';jac_c1,zeros(2,2)];
    Bf = [B_f*u-C_f*dq_f-G_f;-gamma_c1];
    temp = Af\Bf;
    ddq = temp(1:7);
   
    %[cmh,cmv,dcmh,dcmv,acmh,acmv] = center_of_mass_5DoF([q',dq'],ddq(1:5))
    
else % double support
end
dx = [dq_f;ddq];


% --------------------------------------------------------------------------

function [value,isterminal,direction] = events(t,x)
% Locate the time when critial angle of stance leg minus stance leg angle
% passes through zero in a decreasing direction and stop integration.
x_red = [x(1:5);x(8:12)];
if 1
  % when swing leg foot touches ground
  value(1) = swing_foot_height(x_red');
    
  isterminal(1) = 1;  % stop when this event occurs
  direction(1) = -1;  % decreasing direction detection
else
  %% noevents, just watch him fall!!
  value=1;
  isterminal=1;
  direction=1;
end

% --------------------------------------------------------------------------

function [tout,xout_f,torqueout,yout,dyout,fout,step_t] = demo
global optim_var

torque = [];
t_2 = [];
y = [];
dy = [];
force = [];

tstart = 0;
tfinal = 30;

x0 = initialize_HZD(optim_var) ;% Run intialization routine HZD

%x0 = [];
options = odeset('Events','on','Refine',4,'RelTol',10^-5,'AbsTol',10^-6);

tout = tstart;
xout = x0.';
xout_f = xout; % 7 DoF
xout = [xout(:,1:5),xout(:,8:12)];
torqueout = control_synth(xout,optim_var);
teout  = []; xeout = []; ieout = [];
step_t = [];
[h,lfh] = decouple_5links_5DoF(xout(1:5)',xout(6:10)',optim_var); % for controller
yout = [y; h.'];
dyout = [dy; lfh.'];
dxout = f(0,xout_f');
[f_tan,f_norm] = stance_force(xout',dxout(8:12)');
fout = [f_tan f_norm];

modelP; % get model parameters
numSteps = 1;
for i = 1:1:numSteps % steps
  % Solve until the first terminal event.
  [t,x,te,xe,ie] = ode45('walk_HZD',[tstart tfinal],x0,options);

  nt = length(t);
  
  % Accumulate output.  This could be passed out as output arguments.
  tout = [tout; t(2:nt)];
  xout_f = [xout_f; x(2:nt,:)];
  torque = control_synth([x(2:nt,1:5),x(2:nt,8:12)],optim_var);
  torqueout =[torqueout; torque];
  y = [];
  dy = [];
  force = [];
  xr = [x(1:nt,1:5),x(1:nt,8:12)];
  [row,col] = size(xr(2:nt,:));
  
  for k = 1:row
    [h,lfh] = decouple_5links_5DoF(xr(k,1:5)',xr(k,6:10)',optim_var);
    dx = f(0,x(k,:)');
    [f_tan,f_norm] = stance_force(xr(k,:),dx(8:12));
    force = [force; f_tan f_norm];
    y = [y; h.'];
    dy = [dy; lfh.'];
  end
  yout = [yout; y];
  dyout = [dyout; dy];
  fout = [fout; force];

  if t(nt) <= tfinal
     
    % Set the new initial conditions (after impact).
    x0 = impact_5links_5DoF_HZD(x(nt,:));
    
    step_len = step_length(xr(nt,:));

    % Diplay some useful info
    im = sprintf('%.5f',x0(16));
    dt = sprintf('%.5f',t(length(t))-t(1));
    sl = sprintf('%.5f',step_len);
    v = step_len/(t(length(t))-t(1));
    v_txt = sprintf('%.5f',v);
    disp(['step: ',num2str(i),', impact: ',im, ...
	  ', delta_t: ',dt,', length: ',sl,', rate: ',v_txt]);

    step_t = [step_t t(length(t))];
    
    x0 = x0(1:14);
    
    % Stop simulation if swing foot is behind stance foot at end of "step"
    if step_len <= 0
      disp('Stepping error ... swing foot didn''t swing forward.');
    end

    [pHh,pHv] = hip_pos(xr(nt,:));
    
    % Stop simulation if hips are too low
    if pHv <= L3
      disp('Stepping error ... hips too low.');
    end
  end
  
  tstart = t(nt);
  if (tstart >= tfinal)
    break
  end
end
