function varargout = walk_VPP(t,x,flag)


if nargin == 2, flag = ''; end

if nargin == 0, flag = 'demo'; end

switch flag
 case ''                                 % Return dx/dt = dynamics(t,x).
  varargout{1} = f(t,x);
 case 'events'                           % Return [value,isterminal,direction].
  [varargout{1:3}] = events(t,x);
 case 'demo'                             % Run a demo.
  [varargout{1:6}] = demo;
 otherwise
  error(['Unknown flag ''' flag '''.']);
end


% --------------------------------------------------------------------------
% This is the system dynamics function

function dx = f(~,x)
%% Unwrap Parameter
%gamma_vpp = param.gamma_vpp;    % VPP orientation from centerline
%r_vpp = param.r_vpp;            % COM to VPP
%r_h = param.r_h;                % Hip to COM

% Leg stiffness
%q = x(1:5); dq = x(8:12);
q_f = x(1:7);
dq_f = x(8:14);

[D_f,C_f,G_f,B_f] = dyn_mod_5links_7DoF(q_f,dq_f);
[jac_c1,jac_c2,gamma_c1,gamma_c2,~,~] = constraint_5links_7DoF(q_f,dq_f);
     
        % the VPP controller
        % return the geometry of robot
        u = Controller_VPP_store(q_f,dq_f);
       
%ddq = inv_D*(-C*dq-G+B*u);
% D*ddq + C*dq + G = B*u +J'*F
FSM = 1;

if FSM ==1 % Single support
    Af = [D_f,-jac_c1';jac_c1,zeros(2,2)];
    Bf = [B_f*u-C_f*dq_f-G_f;-gamma_c1];
    temp = Af\Bf;
    ddq = temp(1:7);
    %jac_c1*dq_f
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

function [tout,xout_f,torqueout,fout,step_t,VBLAout] = demo
global param 
torque = [];
VBLA = [];
t_2 = [];

force = [];

tstart = 0;
tfinal = 30;

x0 = initialize_VPP() ;% Run intialization routine HZD
%x0 = [];
options = odeset('Events','on','Refine',4,'RelTol',10^-5,'AbsTol',10^-6);

tout = tstart;
xout = x0.';
xout_f = xout; % 7 DoF
xout = [xout(:,1:5),xout(:,8:12)];
[torqueout,~,VBLAout] = Controller_VPP_store(xout_f(:,1:7)',xout_f(:,8:14)');
torqueout = torqueout';
teout  = []; xeout = []; ieout = [];
step_t = [];

dxout = f(0,xout_f');
GRF = stance_force_VPP(xout_f(:,1:7)',xout_f(:,8:14)');
fout = [GRF(1) GRF(2)];

modelP; % get model parameters
numSteps = param.numSteps;

for i = 1:1:numSteps % steps
%% 
    %{
  if ( i == 11)
     param.offset = pi/50; 
  elseif ( i ==15 )
      param.offset = pi/20;
  elseif ( i ==19)
     param.offset = pi/35;
  end
   %}
%%
  % Solve until the first terminal event.
  [t,x,te,xe,ie] = ode45('walk_VPP',[tstart tfinal],x0,options);

  nt = length(t);
  
  % Accumulate output.  This could be passed out as output arguments.
  tout = [tout; t(2:nt)];
  xout_f = [xout_f; x(2:nt,:)];
  
  [~,torque,VBLA] = Controller_VPP_store(xout_f(2:nt,1:7)',xout_f(2:nt,8:14)');
  
  torqueout =[torqueout; torque];
  VBLAout = [VBLAout;VBLA];
  force = [];
  xr = [x(1:nt,1:5),x(1:nt,8:12)];
  
  
  GRF = stance_force_VPP(xout_f(2:nt,1:7)',xout_f(2:nt,8:14)');
 
  fout = [fout; GRF];

  if t(nt) <= tfinal
     
    % Set the new initial conditions (after impact).
    x0 = impact_5links_5DoF_VPP(x(nt,:));
    
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
