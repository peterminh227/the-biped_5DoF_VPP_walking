function [ST_length_T_S,ST_theta,ST_J_polar_T_S,SW_length_T_S,SW_theta,SW_J_polar_T_S] = jacobian_stance(q,dq)
% JACOBIAN_STANCE
%   [ST_LENGTH_T_S,ST_THETA,ST_J_POLAR_T_S,SW_LENGTH_T_S,SW_THETA,SW_J_POLAR_T_S] = JACOBIAN_STANCE(Q,DQ)  mapping the BTSLIP to the rigid body

%Eric Westervelt
%2016 Version: Peter Minh
%06-Dec-2016 00:29:49

modelP;

q31L=q(1); q32L=q(2); q41L=q(3); q42L=q(4); q1L=q(5);
dq31L=dq(1); dq32L=dq(2); dq41L=dq(3); dq42L=dq(4); dq1L= dq(5);


% ST_theta, ST_length_T_S
ST_length_T_S = (abs(L4*sin(q41L))^2 + abs(L3 + L4*cos(q41L))^2)^(1/2);
ST_theta = q31L - pi/2 + atan2(-(L3*cos(q31L + q1L)^2 + L3*sin(q31L + q1L)^2 + L4*cos(q31L + q41L + q1L)*cos(q31L + q1L) + L4*sin(q31L + q41L + q1L)*sin(q31L + q1L))/(cos(q31L + q1L)^2 + sin(q31L + q1L)^2), -(L4*cos(q31L + q41L + q1L)*sin(q31L + q1L) - L4*sin(q31L + q41L + q1L)*cos(q31L + q1L))/(cos(q31L + q1L)^2 + sin(q31L + q1L)^2));

% ST_J_polar_T_S
ST_J_polar_T_S=zeros(2,7);
ST_J_polar_T_S(1,3)=(2*L4*sign(L4*sin(q41L))*cos(q41L)*abs(L4*sin(q41L)) - 2*L4*abs(L3 + L4*cos(q41L))*sign(L3 + L4*cos(q41L))*sin(q41L))/(2*(abs(L4*sin(q41L))^2 + abs(L3 + L4*cos(q41L))^2)^(1/2));
ST_J_polar_T_S(2,1)=1;
ST_J_polar_T_S(2,3)=(2*imag(L4*cos(q41L))*real(L4*sin(q41L)) - 2*real(L4*cos(q41L))*imag(L4*sin(q41L)) + imag(L4*cos(q41L))^2 + real(L4*cos(q41L))^2 + imag(L4*sin(q41L))^2 + real(L4*sin(q41L))^2 + imag(L3)*imag(L4*cos(q41L)) + real(L3)*real(L4*cos(q41L)) + imag(L3)*real(L4*sin(q41L)) - real(L3)*imag(L4*sin(q41L)))/((imag(L3) + imag(L4*cos(q41L)))^2 - 2*imag(L4*sin(q41L))*(real(L3) + real(L4*cos(q41L))) + (real(L3) + real(L4*cos(q41L)))^2 + imag(L4*sin(q41L))^2 + real(L4*sin(q41L))^2 + 2*real(L4*sin(q41L))*(imag(L3) + imag(L4*cos(q41L))));

% SW_theta, SW_length_T_S
SW_length_T_S = (abs(L4*sin(q42L))^2 + abs(L3 + L4*cos(q42L))^2)^(1/2);
SW_theta = q32L - pi/2 + atan2(-(L3*cos(q32L + q1L)^2 + L3*sin(q32L + q1L)^2 + L4*cos(q32L + q42L + q1L)*cos(q32L + q1L) + L4*sin(q32L + q42L + q1L)*sin(q32L + q1L))/(cos(q32L + q1L)^2 + sin(q32L + q1L)^2), -(L4*cos(q32L + q42L + q1L)*sin(q32L + q1L) - L4*sin(q32L + q42L + q1L)*cos(q32L + q1L))/(cos(q32L + q1L)^2 + sin(q32L + q1L)^2));

% SW_J_polar_T_S
SW_J_polar_T_S=zeros(2,7);
SW_J_polar_T_S(1,4)=(2*L4*sign(L4*sin(q42L))*cos(q42L)*abs(L4*sin(q42L)) - 2*L4*abs(L3 + L4*cos(q42L))*sign(L3 + L4*cos(q42L))*sin(q42L))/(2*(abs(L4*sin(q42L))^2 + abs(L3 + L4*cos(q42L))^2)^(1/2));
SW_J_polar_T_S(2,2)=1;
SW_J_polar_T_S(2,4)=(2*imag(L4*cos(q42L))*real(L4*sin(q42L)) - 2*real(L4*cos(q42L))*imag(L4*sin(q42L)) + imag(L4*cos(q42L))^2 + real(L4*cos(q42L))^2 + imag(L4*sin(q42L))^2 + real(L4*sin(q42L))^2 + imag(L3)*imag(L4*cos(q42L)) + real(L3)*real(L4*cos(q42L)) + imag(L3)*real(L4*sin(q42L)) - real(L3)*imag(L4*sin(q42L)))/((imag(L3) + imag(L4*cos(q42L)))^2 - 2*imag(L4*sin(q42L))*(real(L3) + real(L4*cos(q42L))) + (real(L3) + real(L4*cos(q42L)))^2 + imag(L4*sin(q42L))^2 + real(L4*sin(q42L))^2 + 2*real(L4*sin(q42L))*(imag(L3) + imag(L4*cos(q42L))));
