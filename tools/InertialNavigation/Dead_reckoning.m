function [p_DR,v_DR,q_DR,P_DR]=Dead_reckoning(N,T,y_gyr,y_acc,g,q_0,v_0,p_0,R_pos,R_acc,R_gyr)
%Copyright (C) 2022 by Frida Viset

%%Inputs
%N - number of samples
%T - sampling time
%y_gyr - gyroscope measurements, 3 times N
%y_acc - accelerometer measurements, 3 times N
%g - gravity vector, in world frame, 3 times N
%q_0 - initial orientation, unit quaternion, 4 times 1
%v_0 - initial velocity, 3 times 1
%p_0 - initial position, 3 times 1

%%Outputs
%p_DR - position estimate, 3 times N
%v_DR - velocity estimate, 3 times N
%q_DR - orientation estimate, unit quaternion, 4 times N

%Assign memory to states
p_DR=zeros(3,N);
v_DR=zeros(3,N);
q_DR=zeros(4,N);
f_b_DR=zeros(3,N);
acc_DR=zeros(3,N);

%Set initial orientation
q_DR(:,1)=q_0;
v_DR(:,1)=v_0;
p_DR(:,2)=p_0;

%Set the initial covariance
P_DR=zeros(9,9,N);
P_DR(:,:,1)=0.0001*eye(9);

for t=2:N
    %DR update
    f_b_DR(:,t-1)=y_acc(:,t-1);
    acc_DR(:,t-1)=quat2Rot(q_DR(:,t-1))*f_b_DR(:,t-1)+g;
    v_DR(:,t)=v_DR(:,t-1)+T*acc_DR(:,t-1);
    p_DR(:,t)=p_DR(:,t-1)+T*v_DR(:,t-1);
    q_DR(:,t)=exp_q_R(y_gyr(:,t-1)*T,q_DR(:,t-1));
    
    
    %Covariance of DR estimate
    yw=quat2Rot(q_DR(:,t-1))*y_acc(:,t-1);
    ywx=[0, -yw(3), yw(2);
        yw(3), 0, -yw(1);
        -yw(2), yw(1), 0];
    G_t=[eye(3), zeros(3), zeros(3);
        zeros(3), T*quat2Rot(q_DR(:,t)), zeros(3);
        zeros(3), zeros(3), T*quat2Rot(q_DR(:,t))];
    F_t=[eye(3), T*eye(3), zeros(3);
        zeros(3), eye(3), -T*ywx;
        zeros(3), zeros(3), eye(3)];
    P_DR(:,:,t)=F_t*P_DR(:,:,t-1)*F_t'+G_t*[R_pos, zeros(3), zeros(3); zeros(3), R_acc, zeros(3); zeros(3), zeros(3), R_gyr]*G_t';
    P_DR(:,t)=P_DR(:,t-1);
end

end