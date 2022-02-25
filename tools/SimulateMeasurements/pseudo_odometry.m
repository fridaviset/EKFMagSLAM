function [delta_p,delta_q,p_DR,q_DR,N]=pseudo_odometry(p,q,R_p,R_q,p_0,q_0,T,bias)
%Copyright (C) 2022 by Frida Viset

%Create noisy odometry measurements of position
N=size(p,2);
delta_p=[p(:,2:end)-p(:,1:end-1),zeros(3,1)]+mvnrnd(zeros(3,1),R_p,N)'+bias;

%Deduce some orientation odometry
delta_q=zeros(3,N);
for t=1:N-1
    q_t_C=[-q(1,t); q(2:4,t)];
    delta_q(:,t)=quat2angleaxis(quatprod(q(:,t+1),q_t_C))+mvnrnd(zeros(3,1),R_q)';
end

%Deduce the position and orientation DR
p_DR=zeros(3,N);
q_DR=zeros(4,N);
p_DR(:,1)=p_0;
q_DR(:,1)=q_0;

for t=2:N
    p_DR(:,t)=p_DR(:,t-1)+delta_p(:,t-1);
    q_DR(:,t)=exp_q_L(delta_q(:,t-1),q_DR(:,t-1));
end


end