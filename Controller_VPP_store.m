function [u,u_store,VBLA] = Controller_VPP_store(q_f,dq_f)

global param
alpha0 = param.alpha0;          % Angle of Attack at Landing
L0 = param.L0;                  % Rest Length
k = param.k;  
c = param.c;
offset = param.offset;
q_f_store = q_f;
dq_f_store = dq_f;
u_store = [];
VBLA = [];
        for i = 1:size(q_f,2)
        
        q = q_f_store(1:5,i);
        dq = dq_f_store(1:5,i);
        %q_f = q_f_store(i,:);
        dq_f = dq_f_store(:,i)';
        out = limb_position(q',0); 
        % change co-ordinate to vector
        [ST_length_T_S,ST_theta,ST_J_polar_T_S,SW_length_T_S,SW_theta,SW_J_polar_T_S] = jacobian_stance(q,dq);
        
        % stance leg 
        ST_L = ST_length_T_S;
        ST_J = ST_J_polar_T_S;
        ST_dX = ST_J*dq_f';
        ST_dL = ST_dX(1);
        ST_dth = ST_dX(2);
        % swing leg 
        SW_L = SW_length_T_S;
        SW_J = SW_J_polar_T_S;
        SW_dX = SW_J*dq_f';
        SW_dL = SW_dX(1);
        SW_dth = SW_dX(2);
        % stance leg control
        ST_Fs = k*(L0-ST_L) -k/100*ST_dL;
        if ST_Fs < 0 
            ST_Fs = 0 - k/100*ST_dL;
        end
        
        % Ft by JW innovation
        % define vector leg
        ST_B_L_C = [out.pCoM1 - out.pFoot11;out.pCoM2 - out.pFoot12];
        ST_B_L_H = [out.pH1 - out.pFoot11;out.pH2 - out.pFoot12];
        ST_zeta = angle2vec(ST_B_L_C,ST_B_L_H);
        cur_phi = q(5)-offset; % trunk angle
        ST_cur_beta = -c*cur_phi; % ST_cur_beta*cur_phi >= 0
        ST_beta = ST_zeta + ST_cur_beta;
        ST_Ft = ST_Fs*tan(ST_beta);
        ST_Tau_ff = ST_L*ST_Ft;
        ST_F = [ST_Fs;-ST_Tau_ff];
        %% VBLA calculation for swing leg angle
        V = [dq_f(6);dq_f(7)]*1/sqrt(9.81*ST_length_T_S);
        G = [0;-1];
        muy = 0.4;
        O = muy*V+(1-muy)*G;
        angle_VBLA = pi + atan2(O(2),O(1));
        VBLA = [VBLA;angle_VBLA];
        %%
        % swing leg control
        enable = 0;
        if ((swing_foot_height(q') < 1e-10) && (enable == 1))
        %% control as the stance phase
        SW_Fs = k*(L0-SW_L) -k/100*SW_dL;
        if SW_Fs < 0 
            SW_Fs = 0 - k/100*SW_dL;
        end    
        SW_B_L_C = [out.pCoM1 - out.pFoot21;out.pCoM2 - out.pFoot22];
        SW_B_L_H = [out.pH1 - out.pFoot21;out.pH2 - out.pFoot22];
        SW_zeta = angle2vec(SW_B_L_C,SW_B_L_H);
        cur_phi = q(5)+offset; % trunk angle
        SW_cur_beta = -c*cur_phi; % ST_cur_beta*cur_phi >= 0
        SW_beta = SW_zeta + SW_cur_beta;
        SW_Ft = SW_Fs*tan(SW_beta);
        SW_Tau_ff = SW_L*SW_Ft;
        SW_F = [SW_Fs;-SW_Tau_ff];
        else
        phi = wrapToPi(cur_phi) + pi/2;
        SW_alpha =  wrapToPi(phi + SW_theta);   
        if ((SW_alpha - pi/2) < deg2rad(10))
            SW_L_d = L0*4/5;
        else
            SW_L_d = L0;
        end
        %SW_F = [0;0];
        SW_Fs = k*(SW_L_d-SW_L) - 100*SW_dL;
        SW_Tau_ff = 0;
        SW_Tau_fb = -5*2*(alpha0 - SW_alpha) + 1*SW_dth;
        SW_F = [SW_Fs;-(SW_Tau_ff+SW_Tau_fb)];
        end
        % torque mapping
        ST_tor = ST_J'*ST_F;
        SW_tor = SW_J'*SW_F;
        D_fric = -diag([0.1 0.1 0.1 0.1]);
        u = [ST_tor(1);SW_tor(2);ST_tor(3);SW_tor(4)] +  D_fric*dq(1:4);
        u_store = [u_store;u'];
        end
        
end