function GRF_stance = stance_force_VPP(q_f,dq_f)
q_f_store = q_f;
dq_f_store = dq_f;
GRF_stance = [];
    for i = 1:size(q_f,2)
         %q = q_f_store(1:5,i);
         %dq = dq_f_store(1:5,i);
         q_f = q_f_store(:,i);
         dq_f = dq_f_store(:,i);
         [D_f,C_f,G_f,B_f] = dyn_mod_5links_7DoF(q_f,dq_f);
         [jac_c1,jac_c2,gamma_c1,gamma_c2,~,~] = constraint_5links_7DoF(q_f,dq_f);
         [u,~] = Controller_VPP_store(q_f,dq_f);
         Af = [D_f,-jac_c1';jac_c1,zeros(2,2)];
         Bf = [B_f*u-C_f*dq_f-G_f;-gamma_c1];
         temp = Af\Bf;
         GRF_stance = [GRF_stance;temp(8) temp(9)];
    end
end