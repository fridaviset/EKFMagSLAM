function [p_EKF,q_EKF,m_EKF,P_EKF,P_pose_prior,P_pose_posterior,runtime]=EKF_quick(N,delta_p,delta_q,y_mag,q_0,p_0,R_p,R_q,sigma_y,Lambda,Indices,sigma_lin,N_m,xl,xu,yl,yu,zl,zu)
%Copyright (C) 2022 by Frida Viset

%Start timer
tic;

%Pre-allocate position trajectory for all particles
q_EKF=zeros(4,N);
p_EKF=zeros(3,N);

%Pre-allocate space for the position covariance
P_EKF=zeros(N_m+9,N_m+9);
%P_EKF(4:6,4:6)=0.001*eye(3);

%Pre-allocate storage for pose covariance estimate
P_pose_prior=zeros(6,6,N);
P_pose_posterior=zeros(6,6,N);

%Initialise orientation and position estimates
q_EKF(:,1)=q_0;
p_EKF(:,1)=p_0;

%Initialise the magnetic field map
P_EKF(7:end,7:end)=[sigma_lin^2*eye(3), zeros(3,N_m);
    zeros(N_m,3), Lambda];
m_EKF=zeros(N_m+3,1);

%iterate through all timesteps
for t=2:N
    
    %KF Dyn update
    p_EKF(:,t)=p_EKF(:,t-1)+delta_p(:,t-1);
    q_EKF(:,t)=exp_q_L(delta_q(:,t-1),q_EKF(:,t-1));
    P_EKF(1:6,1:6)=P_EKF(1:6,1:6)+[R_p, zeros(3); zeros(3), R_q];
    
    %Store pose prior covariances
    P_pose_prior(:,:,t)=P_EKF(1:6,1:6);
    
    %KF meas update
    NablaPhi=[eye(3),Nabla_Phi3D(p_EKF(:,t),N_m,xl,xu,yl,yu,zl,zu,Indices)];
    f=NablaPhi*m_EKF;
    J=reshape(JacobianPhi3D(p_EKF(:,t),N_m,xl,xu,yl,yu,zl,zu,Indices),9,N_m)*m_EKF(4:end);
    J = reshape(J,3,3);
    
    %Prepare scew-symmetric magnetic field vector
    fx=[0 -f(3) f(2);
        f(3) 0 -f(1);
        -f(2) f(1) 0];
    
    %Meas state update full EKF
    H_t=[J, fx, NablaPhi];
    S_t=H_t*P_EKF*H_t'+sigma_y^2*eye(3);
    K_t=P_EKF*H_t'*inv(S_t);
    S_t=0.5*(S_t+S_t');
    eps_t=quat2Rot(q_EKF(:,t))*y_mag(:,t)-(f);
    eta_t=K_t*eps_t;
    p_EKF(:,t)=p_EKF(:,t)+eta_t(1:3);
    q_EKF(:,t)=exp_q_L(eta_t(4:6),q_EKF(:,t));
    m_EKF=m_EKF+eta_t(7:end);
    
    P_EKF=P_EKF-K_t*S_t*K_t';
    P_EKF=1/2*(P_EKF+P_EKF');
    
    %Store pose posterior covariances
    P_pose_posterior(:,:,t)=P_EKF(1:6,1:6);
    
end

runtime=toc;

end