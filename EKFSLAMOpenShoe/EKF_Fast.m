function [p_hat,q_hat,m,P_m,P_pose_prior,P_pose_posterior]=EKF(N,T,foldername,delta_p,delta_q,y_mag,q_0,p_0,R_p,R_q,sigma_y,Lambda,Indices,sigma_lin,N_m,xl,xu,yl,yu,zl,zu,fontsize)
%Copyright (C) 2022 by Frida Viset

%Pre-allocate position trajectory for all particles
q_hat=zeros(4,N);
p_hat=zeros(3,N);

%Pre-allocate storage for pose covariance estimate
P_pose_prior=zeros(6,6,N);
P_pose_posterior=zeros(6,6,N);

%Pre-allocate space for the position covariance
P=zeros(N_m+9,N_m+9);
P(1:6,1:6)=0.001*eye(6);

%Initialise orientation and position estimates
q_hat(:,1)=q_0;
p_hat(:,1)=p_0;

%Use the limits of the magnetic field norm to have reasonable limits for
%the plot
y_mag_norm=sqrt(sum(y_mag.^2,1));
cmin=min(y_mag_norm);
cmax=max(y_mag_norm);

%Initialise the magnetic field map
P(7:end,7:end)=[sigma_lin^2*eye(3), zeros(3,N_m);
    zeros(N_m,3), Lambda];
m=zeros(N_m+3,1);

%iterate through all timesteps
for t=2:N
    
    %KF Dyn update
    p_hat(:,t)=p_hat(:,t-1)+delta_p(:,t-1);
    q_hat(:,t)=exp_q_L(delta_q(:,t-1),q_hat(:,t-1));
    
    P(1:6,1:6)=P(1:6,1:6)+[R_p, zeros(3); zeros(3), R_q];
    
    %Store pose prior covariances
    P_pose_prior(:,:,t)=P(1:6,1:6);
    
    if (isnan(y_mag(1,t)))
        %Skip meas update
    else
        
        %KF meas update
        NablaPhi=[eye(3),Nabla_Phi3D(p_hat(:,t),N_m,xl,xu,yl,yu,zl,zu,Indices)];
        f=NablaPhi*m;
        J=reshape(JacobianPhi3D(p_hat(:,t),N_m,xl,xu,yl,yu,zl,zu,Indices),9,N_m)*m(4:end);
        J = reshape(J,3,3);
        
        %Prepare scew-symmetric magnetic field vector
        fx=[0 -f(3) f(2);
            f(3) 0 -f(1);
            -f(2) f(1) 0];
        
        %Meas state update full EKF
        H_t=[J, fx, NablaPhi];
        S_t=H_t*P*H_t'+sigma_y^2*eye(3);
        K_t=P*H_t'*inv(S_t);
        S_t=0.5*(S_t+S_t');
        eps_t=quat2Rot(q_hat(:,t))*y_mag(:,t)-(f);
        eta_t=K_t*eps_t;
        p_hat(:,t)=p_hat(:,t)+eta_t(1:3);
        q_hat(:,t)=exp_q_L(eta_t(4:6),q_hat(:,t));
        m=m+eta_t(7:end);
        
        P=P-K_t*S_t*K_t';
        P=1/2*(P+P');
        
        %Store pose posterior covariances
        P_pose_posterior(:,:,t)=P(1:6,1:6);
    end

    
    P_m=P(7:end,7:end);

end