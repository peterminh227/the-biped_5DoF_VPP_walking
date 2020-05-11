function varargout = sim_zd(t,x,flag)
%SIM_ZD   Simulate the zero dynamics of the kneed biped walker.

%Copyright (c) 2003 by Eric Westervelt and Jessy Grizzle.  This code
%may be freely used for noncommercial ends. If use of this code in
%part or in whole results in publication, proper citation must be
%included in that publication.  This code comes with no guarantees or
%support.

%Eric R. Westervelt & Jessy W. Grizzle
%10/26/00
%11/30/03 - updated for the CDC

if nargin == 2, flag = ''; end

if nargin == 0, flag = 'demo'; end
switch flag
 case ''                                 % Return dx/dt = dynamics(t,x).
  varargout{1} = f(t,x);
 case 'events'                           % Return [value,isterminal,direction].
  [varargout{1:3}] = events(t,x);
 case 'demo'                             % Run a demo.
  [varargout{1:8}] = demo;
 otherwise
  error(['Unknown flag ''' flag '''.']);
end



% --------------------------------------------------------------------------
% This is the system dynamics function

function dz = f(t,z)
global optim_var
dz = zd(z,optim_var)';

% --------------------------------------------------------------------------
% This is the events function

function [value,isterminal,direction] = events(t,z)
%Locate the time when critial angle of stance leg minus stance leg
%angle passes through zero in a decreasing direction and stop
%integration.
global optim_var
x = zd_s_lift(z,optim_var);

[pHh,pHv]=hip_pos(x);

h = swing_foot_height(x');

value(1) = h;
value(2) = pHv - .6;   % hips not too low
value(3) = x(3) - pi;  % stance knee not too bent
value(4) = x(4) - pi; % swing knee not too bent

isterminal(1) = 1;   % stop when this event occurs
isterminal(2) = 1;   % stop  "
isterminal(3) = 1;   % stop  "
isterminal(4) = 1;   % stop  "

direction(1) = -1;  % decreasing direction
direction(2) = 0;   % decreasing  "
direction(3) = 0;   % increasing  "
direction(4) = 0;   % decreasing  "

% --------------------------------------------------------------------------
% This is the demo function

function [tout,zout,step_t,teout,zeout,ieout,J,step_len] = demo
global optim_var

x0 = initialize(optim_var); %

% project initial condition onto zero dynamics
z0 = zd_project(x0);

num_steps = 1;
zout = z0;

tstart = 0;
tfinal = 2;

tout = tstart;

step_t = [];

modelP; % get model parameters

options = odeset('Events','on', ...
		 'Refine',4, ...
		 'RelTol',10^-6, ...
		 'AbsTol',10^-7);

teout = [];
zeout = [];
ieout = [];

for i = 1:num_steps

  % Solve until the first terminal event.
  [t,z,te,ze,ie] = ode45('sim_zd',[tstart tfinal],z0,options);
  nt = length(t);
  
  % Accumulate output.  This could be passed out as output
  % arguments.
  tout = [tout; t(2:nt)];
  zout = [zout; z(2:nt,:)];
  teout = [teout; te];
  zeout = [zeout; ze];
 ieout = [ieout; ie];
  if t(nt) <= tfinal %& i ~= num_steps
    % Set the new initial conditions (after impact).
    xf = zd_s_lift(z(nt,:),optim_var);
    x0 = impact_5links_5DoF(xf);
    z0 = zd_project(x0(1:10)');

    step_len = step_length(xf);
    
    % Diplay some useful info
    im = sprintf('%.5f',x0(12));
    dt = sprintf('%.5f',t(length(t))-t(1));
    sl = sprintf('%.5f',step_len);
    v = step_len/(t(length(t))-t(1));
    v_txt = sprintf('%.5f',v);
    J = 1/(step_len)*z(nt,3);
    %J = 1/(step_len*9.81*(12+6.8+3.2))*z(nt,3);
    J_txt = sprintf('%.5f',J);
    disp(['step: ',num2str(i),', impact: ',im, ...
	  ', delta_t: ',dt,', length: ',sl,', rate: ',v_txt,', Cost of trans: ',J_txt]);
    
    step_t = [step_t t(length(t))];

    % Stop simulation if swing foot is behind stance foot at end of
    % "step"
    if step_len <= 0
     % disp('stepping error - swing foot didn''t swing forward.');
    end
  end

  tstart = t(nt);
  STEP_TIME = t(nt);
  if tstart >= tfinal
    break
  end  
end

