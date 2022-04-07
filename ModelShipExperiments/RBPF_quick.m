function [p_RBPF,q_RBPF,runtime]=RBPF_quick(N,delta_p,delta_q,y_mag,q_0,p_0,R_p,R_q,sigma_y,Lambda,Indices,sigma_lin,N_m,xl,xu,yl,yu,zl,zu,N_p)
%Copyright (C) 2022 by Frida Viset

%Start timer
tic;

%Pre-allocate position trajectory for all particles
q_cloud=zeros(4,N,N_p);
p_cloud=zeros(3,N,N_p);

%Pre-allocate the max likelihood estimates
p_RBPF=zeros(3,N);
p_RBPF(:,1)=p_0;

%Pre-allocate variances of the weighted mean estimates
P_hat=zeros(3,3,N);

%Initialise orientation and position estimates
for i=1:N_p
    q_cloud(:,1,i)=q_0;
    p_cloud(:,1,i)=p_0;
end

%Initialise weights
weights=ones(N_p,1)./N_p;

%Initialise the magnetic field map
P_m=zeros(N_m+3,N_m+3,N_p);
m=zeros(N_m+3,N_p);
for i=1:N_p
    P_m(:,:,i)=[sigma_lin^2*eye(3), zeros(3,N_m);
        zeros(N_m,3), Lambda];
end

%Allocate space for storing estimated measurements before and after resampling
y_prior=zeros(3,N_p);
y_post=zeros(3,N_p);

%Iterate through all timesteps
for t=2:N
    
    %PF Dyn update
    for i=1:N_p
        e_p_sim=mvnrnd(zeros(3,1),R_p)';
        e_q_sim=mvnrnd(zeros(3,1),R_q)';
        p_cloud(:,t,i)=p_cloud(:,t-1,i)+delta_p(:,t-1)+e_p_sim;
        q_cloud(:,t,i)=exp_q_L(delta_q(:,t-1)+e_q_sim,q_cloud(:,t-1,i));
    end
    
    %KF Dyn update - none, because the magnetic field is assumed stationary
    
    %PF Meas update
    for i=1:N_p
        NablaPhi=[eye(3),Nabla_Phi3D(p_cloud(:,t,i),N_m,xl,xu,yl,yu,zl,zu,Indices)];
        f=NablaPhi*m(:,i);
        H_t=NablaPhi;
        S_t=H_t*P_m(:,:,i)*H_t'+sigma_y^2*eye(3);
        S_t=0.5*(S_t+S_t');
        weights(i)=mvnpdf(quat2Rot(q_cloud(:,t,i))*y_mag(:,t),f,S_t)*weights(i);
        y_prior(:,i)=quat2Rot(q_cloud(:,t,i))'*f;
    end
    
    [q_cloud,p_cloud]=resample(q_cloud,p_cloud,weights);
    weights=0*weights+1./N_p;
    
    
    %KF Meas update
    for i=1:N_p
        %KF meas update
        NablaPhi=[eye(3),Nabla_Phi3D(p_cloud(:,t,i),N_m,xl,xu,yl,yu,zl,zu,Indices)];
        f=NablaPhi*m(:,i);
        H_t=NablaPhi;
        S_t=H_t*P_m(:,:,i)*H_t'+sigma_y^2*eye(3);
        K_t=P_m(:,:,i)*H_t'*inv(S_t);
        S_t=0.5*(S_t+S_t');
        eps_t=quat2Rot(q_cloud(:,t,i))*y_mag(:,t)-(f);
        eta_t=K_t*eps_t;
        m(:,i)=m(:,i)+eta_t;
        P_m(:,:,i)=P_m(:,:,i)-K_t*S_t*K_t';
        P_m(:,:,i)=1/2*(P_m(:,:,i)+P_m(:,:,i)');
        y_post(:,i)=quat2Rot(q_cloud(:,t,i))'*f;
    end
    
end

[~,arg_max]=max(weights);
p_RBPF=p_cloud(:,:,arg_max);
q_RBPF=q_cloud(:,:,arg_max);

runtime=toc;

end