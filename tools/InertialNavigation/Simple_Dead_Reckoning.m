function [p_DR,q_DR]=Simple_Dead_Reckoning(delta_p,delta_q,p_0,q_0)
%Copyright (C) 2022 by Frida Viset

=size(delta_p,2);
p_DR=zeros(3,N);
q_DR=zeros(4,N);
p_DR(:,1)=p_0;
q_DR(:,1)=q_0;
for t=2:N
    p_DR(:,t)=p_DR(:,t-1)+delta_p(:,t-1);
    q_DR(:,t)=exp_q_R(delta_q(:,t-1),q_DR(:,t-1));
end
end