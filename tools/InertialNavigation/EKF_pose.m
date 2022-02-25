function [p_hat,v_hat,q_hat]=EKF_pose(N,T,y_gyr,y_acc,y_mag,g,q_0,v_0,p_0,R_gyr,R_acc,R_pos,R_mag,Indices,m_sim,N_m,xl,xu,yl,yu,zl,zu,P_0_EKF)
%Copyright (C) 2022 by Frida Viset

%%Inputs
%N - number of samples
%T - sampling time
%y_gyr - gyroscope measurements, 3 times N
%y_acc - accelerometer measurements, 3 times N
%g - gravity field, in world frame, 3 times 1
%m - magnetic field, in world frame, 3 times 1
%q_0 - initial orientation, unit quaternion, 4 times 1
%v_0 - initial velocity, 3 times 1
%p_0 - initial position, 3 times 1
%R_gyr - variance of gyroscope measurements, 3 times 3
%R_acc - variance of accelerometer measurement, 3 times 3

%%Outputs
%p_DR - position estimate, 3 times N
%v_DR - velocity estimate, 3 times N
%q_DR - orientation estimate, unit quaternion, 4 times N

%Assign memory to states
p_hat=zeros(3,N);
v_hat=zeros(3,N);
q_hat=zeros(4,N);
f_b_hat=zeros(3,N);
acc_hat=zeros(3,N);

%Assign memory to covar
P=zeros(9,9,N);

%Initial state deviaton covar
P(:,:,1)=P_0_EKF;

%Prepare scew-symmetric gravity vector
gx=[0 -g(3) g(2);
    g(3) 0 -g(1);
    -g(2) g(1) 0];

%Set initial orientation
p_hat(:,1)=p_0;
v_hat(:,1)=v_0;
q_hat(:,1)=q_0;

%Run the filters
for t=2:N
    
    %Dyn state update full EKF
    f_b_hat(:,t-1)=y_acc(:,t-1);
    acc_hat(:,t-1)=quat2Rot(q_hat(:,t-1))*f_b_hat(:,t-1)+g;
    v_hat(:,t)=v_hat(:,t-1)+T*acc_hat(:,t-1);
    p_hat(:,t)=p_hat(:,t-1)+T*v_hat(:,t-1);
    q_hat(:,t)=exp_q_R(y_gyr(:,t-1)*T,q_hat(:,t-1));
    
    %Dyn covar update full EKF
    yw=quat2Rot(q_hat(:,t-1))*y_acc(:,t-1);
    ywx=[0, -yw(3), yw(2);
        yw(3), 0, -yw(1);
        -yw(2), yw(1), 0];
    G_t=[eye(3), zeros(3), zeros(3);
        zeros(3), T*quat2Rot(q_hat(:,t)), zeros(3);
        zeros(3), zeros(3), T*quat2Rot(q_hat(:,t))];
    F_t=[eye(3), T*eye(3), zeros(3);
        zeros(3), eye(3), -T*ywx;
        zeros(3), zeros(3), eye(3)];
    P(:,:,t)=F_t*P(:,:,t-1)*F_t'+G_t*[R_pos, zeros(3), zeros(3); zeros(3), R_acc, zeros(3); zeros(3), zeros(3), R_gyr]*G_t';
    
    %Acceleration == gravity field assumption
%     if abs(norm(y_acc(:,t))-norm(g))<0.01
%         
%         %Meas state update full EKF
%         H_t=-[zeros(3), zeros(3), (quat2Rot(q_hat(:,t)))'*gx];
%         S_t=H_t*P(:,:,t)*H_t'+R_acc;
%         K_t=P(:,:,t)*H_t'*inv(S_t);
%         eps_t=y_acc(:,t)+quat2Rot(q_hat(:,t))'*(g);
%         eta=K_t*eps_t;
%         q_hat(:,t)=exp_q_L(eta(7:9),q_hat(:,t));
%         v_hat(:,t)=v_hat(:,t)+eta(4:6);
%         p_hat(:,t)=p_hat(:,t)+eta(1:3);
%         
%         %Meas covar update full EKF
%         P(:,:,t)=P(:,:,t)-K_t*S_t*K_t';
%         
%     end
    
    m=Nabla_Phi3D(p_hat(:,t),N_m,xl,xu,yl,yu,zl,zu,Indices)*m_sim;
    J=reshape(JacobianPhi3D(p_hat(:,t),N_m,xl,xu,yl,yu,zl,zu,Indices),9,200)*m_sim;
    J = reshape(J,3,3);
    %Prepare scew-symmetric magnetic field vector
    mx=[0 -m(3) m(2);
        m(3) 0 -m(1);
        -m(2) m(1) 0];
    
    H_t=[J, zeros(3), (quat2Rot(q_hat(:,t)))'*mx];
    S_t=H_t*P(:,:,t)*H_t'+R_mag;
    K_t=P(:,:,t)*H_t'*inv(S_t);
    eps_t=y_mag(:,t)-quat2Rot(q_hat(:,t))'*(m);
    eta=K_t*eps_t;
    q_hat(:,t)=exp_q_L(eta(7:9),q_hat(:,t));
    v_hat(:,t)=v_hat(:,t)+eta(4:6);
    p_hat(:,t)=p_hat(:,t)+eta(1:3);
    
    %Meas covar update full EKF
    P(:,:,t)=P(:,:,t)-K_t*S_t*K_t';
end

end