function u = control_synth(x,optim_var)
%CONTROL_SYNTH    Synthesize control associated with the supplied
%   state.
%
%   U = CONTROL_SYNTH(X)

%Copyright (c) 2003 by Eric Westervelt and Jessy Grizzle.  This code
%may be freely used for noncommercial ends. If use of this code in
%part or in whole results in publication, proper citation must be
%included in that publication.  This code comes with no guarantees or
%support.

%Eric R. Westervelt & Jessy W. Grizzle
%3/13/01
%11/30/03 - updated for the CDC
%>edit by Peter Minh
[n,m] = size(x);

if ~(n == 10 && m == 1) && m ~= 10
  error('input has wrong dimensions');
end

gains = [1, 1, 1, 1];

if ~(n == 10 && m == 1)
  u = zeros(n,4);
  for k = 1:n
    [H,LfH,L2fH,LgLfH] = decouple_5links_5DoF(x(k,1:5)',x(k,6:10)',optim_var);
    
    out = [bb([H(1),LfH(1)],gains(1));
	   bb([H(2),LfH(2)],gains(2));
	   bb([H(3),LfH(3)],gains(3));
	   bb([H(4),LfH(4)],gains(4))];
    
    % Used for controller that use feedback linearization
    u(k,:) = (inv(LgLfH)*(out(:,1)-L2fH))';
  end 
else
  [H,LfH,L2fH,LgLfH] = decouple_5links_5DoF(x(1:5),x(6:10),optim_var);
  
  out = [bb([H(1),LfH(1)],gains(1));
	 bb([H(2),LfH(2)],gains(2));
	 bb([H(3),LfH(3)],gains(3));
	 bb([H(4),LfH(4)],gains(4))];
  
  % Used for controller that use feedback linearization
  u = (inv(LgLfH)*(out(:,1)-L2fH))';
end