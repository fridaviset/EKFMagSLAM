function q_new=exp_q_R(eta,q_old)
%Copyright (C) 2022 by Frida Viset

eta1=eta(1);
eta2=eta(2);
eta3=eta(3);

OMEGA=[0 -eta1 -eta2 -eta3;
    eta1 0 eta3 -eta2;
    eta2 -eta3 0 eta1;
    eta3 eta2 -eta1 0];

normeta=norm(eta);
if normeta~=0
    q_new=(cos(normeta/2)*eye(4)+1/normeta*sin(normeta/2)*OMEGA )*q_old;
    q_new=q_new./norm(q_new); %Normalization for numeric stability
else
    q_new=q_old;
end
end