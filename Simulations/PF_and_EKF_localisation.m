function [p_PF,q_PF,p_EKF,q_EKF]=PF_and_EKF_localisation(N,T,foldername,delta_p,delta_q,y_mag,q_0,p_0,P_0,p,R_p,R_q,m_sim,sigma_y,Indices,N_m,xl,xu,yl,yu,zl,zu,M,fontsize,experiment,param)
%Copyright (C) 2022 by Frida Viset

%Pre-allocate position trajectory for all particles
q_cloud=zeros(4,N,M);
p_cloud=zeros(3,N,M);

%Pre-allocate the weighted mean estimates
p_PF=zeros(3,N);
p_PF(:,1)=p_0;

%Pre-allocate position trajectory for all particles
q_EKF=zeros(4,N);
p_EKF=zeros(3,N);

%Pre-allocate space for the position covariance
P_EKF=P_0;

%Initialise orientation and position estimates
q_EKF(:,1)=q_0;
p_EKF(:,1)=p_0;

%Pre-allocate variances of the weighted mean estimates
P_RBPF=zeros(3,3,N);

%Initialise orientation and position estimates
for i=1:M
    q_cloud(:,1,i)=q_0;
    p_cloud(:,1,i)=p_0+mvnrnd(zeros(3,1),P_0(1:3,1:3))';
end

%Initialise weights
weights=ones(M,1)./M;

%Allocate space for storing estimated measurements before and after resampling
y_prior=zeros(3,M);
y_post=zeros(3,M);

%Iterate through all timesteps
for t=2:N
    
    %KF Dyn update
    p_EKF(:,t)=p_EKF(:,t-1)+delta_p(:,t-1);
    q_EKF(:,t)=exp_q_L(delta_q(:,t-1),q_EKF(:,t-1));
    
    P_EKF(1:6,1:6)=P_EKF(1:6,1:6)+[R_p, zeros(3); zeros(3), R_q];
    
    %KF meas update
    NablaPhi=[eye(3),Nabla_Phi3D(p_EKF(:,t),N_m,xl,xu,yl,yu,zl,zu,Indices)];
    f=NablaPhi*m_sim;
    J=reshape(JacobianPhi3D(p_EKF(:,t),N_m,xl,xu,yl,yu,zl,zu,Indices),9,N_m)*m_sim(4:end);
    J = reshape(J,3,3);
    
    %Prepare scew-symmetric magnetic field vector
    fx=[0 -f(3) f(2);
        f(3) 0 -f(1);
        -f(2) f(1) 0];
    
    %Meas state update full EKF
    H_t=[J, fx];
    S_t=H_t*P_EKF*H_t'+sigma_y^2*eye(3);
    K_t=P_EKF*H_t'*inv(S_t);
    S_t=0.5*(S_t+S_t');
    eps_t=quat2Rot(q_EKF(:,t))*y_mag(:,t)-(f);
    eta_t=K_t*eps_t;
    p_EKF(:,t)=p_EKF(:,t)+eta_t(1:3);
    q_EKF(:,t)=exp_q_L(eta_t(4:6),q_EKF(:,t));
    
    P_EKF=P_EKF-K_t*S_t*K_t';
    P_EKF=1/2*(P_EKF+P_EKF');
    
    %PF Dyn update
    for i=1:M
        e_p_sim=mvnrnd(zeros(3,1),R_p)';
        e_q_sim=mvnrnd(zeros(3,1),R_q)';
        p_cloud(:,t,i)=p_cloud(:,t-1,i)+delta_p(:,t-1)+e_p_sim;
        q_cloud(:,t,i)=exp_q_L(delta_q(:,t-1)+e_q_sim,q_cloud(:,t-1,i));
    end
    
    %KF Dyn update - none, because the magnetic field is assumed stationary
    
    %PF Meas update
    
    for i=1:M
        NablaPhi=[eye(3),Nabla_Phi3D(p_cloud(:,t,i),N_m,xl,xu,yl,yu,zl,zu,Indices)];
        f=NablaPhi*m_sim;
        S_t=sigma_y^2*eye(3);
        S_t=0.5*(S_t+S_t');
        weights(i)=mvnpdf(quat2Rot(q_cloud(:,t,i))*y_mag(:,t),f,S_t)*weights(i);
        y_prior(:,i)=quat2Rot(q_cloud(:,t,i))'*f;
    end
    
    weights=weights./sum(weights);
    [q_cloud,p_cloud]=resample(q_cloud,p_cloud,weights);
    weights=0*weights+1./M;
    
    %Plotting results
    if (t==2 && experiment==1)
        
        h1=figure(1); clf;
        set(h1,'Position',[100 100 500 500])
        [~,arg_max]=max(weights);
        scatter3(p_cloud(1,t,:),p_cloud(2,t,:),p_cloud(3,t,:)+50,max(1,2*M*weights),'ko','linewidth',1);
        hold on;
        plot3(p_EKF(1,t),p_EKF(2,t),p_EKF(3,t)+50,'k+','linewidth',1.2);
        error_ellipse(P_EKF(1:2,1:2),p_EKF(1:2,t));
        plot3(p(1,t),p(2,t),p(3,t)+50,'r+','linewidth',1.2)
        view(2);
        xlabel('x(m)','Interpreter','Latex');
        ylabel('y(m)','Interpreter','Latex');
        res = 0.01;
        z=0;
        [X,Y,Z,pointsVec,PhiVec,NablaPhi3D]=prepare_magnetic_field_plots(xl,xu,yl,yu,zl,zu,N_m,Indices,res,z);
        N_points=size(NablaPhi3D,1);
        NablaPhi3DVec=reshape(NablaPhi3D,N_points*3,N_m+3);
        FieldVec=NablaPhi3DVec*m_sim;
        FieldVec3=reshape(FieldVec,N_points,3);
        N_points=size(FieldVec3,1);
        Norm=zeros(N_points,1);
        for i=1:N_points
            Norm(i)=norm(FieldVec3(i,:));
        end
        Norm=reshape(Norm,size(X));
        surf(X,Y,Z,Norm,'EdgeColor','none');
        view(2);
        xlabel('$x$(m)','Interpreter','Latex','Fontsize',fontsize);
        ylabel('$y$(m)','Interpreter','Latex','Fontsize',fontsize);
        zlabel('$z$(m)','Interpreter','Latex','Fontsize',fontsize);
        c=colorbar;
        axis equal;
        xlim([xl xu]);
        ylim([yl yu]);
        legend({'Particle cloud','EKF estimate','EKF covariance','Ground truth'},'Interpreter','Latex','Fontsize',fontsize);
        
        saveas(gcf,[foldername,'/InitialDistributions/param',num2str(param),'/experiment',num2str(experiment)]);
        saveas(gcf,[foldername,'/InitialDistributions/EPSs/param',num2str(param),'/experiment',num2str(experiment)],'epsc');
    end
    
    %Update the mean state estimate
    x_post=reshape(p_cloud(:,t,:),3,M);
    x_prev=reshape(p_cloud(:,t-1,:),3,M);
    p_PF(:,t)=mean(x_post');
    P_RBPF(:,:,t)=cov(x_post');
    
    %Check the particle deviation
    f=x_prev+delta_p(:,t-1);
    
end

[~,arg_max]=max(weights);
p_PF=p_cloud(:,:,arg_max);
q_PF=q_cloud(:,:,arg_max);


end