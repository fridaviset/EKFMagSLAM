%Copyright (C) 2022 by Frida Viset

close all; clear;

disp('Running Algorithm 1 on all four datasets from model ship');

seed=42;
rng(seed);

experiments=30;
mean_RMSEs_EKF=zeros(experiments,4);
mean_RMSEs_DR=zeros(experiments,4);

for dataset=1:4

for experiment=1:experiments

load(['Aligned_SquareLoop_',num2str(dataset),'.mat']);

R_p=0.0000001*eye(3);
sigma_xy=0.00001;
R_p(1:2,1:2)=sigma_xy*eye(2);
R_q=0.00001*eye(3);
bias=[0.003; 0.003; 0];

%Find room
margin=1;
xl=min(p(1,:))-margin;
xu=max(p(1,:))+margin;
yl=min(p(2,:))-margin;
yu=max(p(2,:))+margin;
zl=min(p(3,:))-margin;
zu=max(p(3,:))+margin;

%Number of basis functions used in Reduced-Rank approximation
N_m=500;

%Magnetic field params
sigma_SE=1;
l_SE=0.8;
sigma_lin=1;
sigma_y=0.1;

%Calculate Lambda and the order of indices used in the
%analytic basis functions of the Reduced-Rank Approximation
[Indices, Lambda]=Lambda3D(N_m,xl,xu,yl,yu,zl,zu,sigma_SE,l_SE);

[delta_p,delta_q,p_DR,q_DR,N]=pseudo_odometry(p,q,R_p,R_q,p_0,q_0,T,bias);

%Adjust odometry noise to compensate for bias
R_p(1:2,1:2)=(0.001+sigma_xy)*eye(2);

%Run the filter
[p_EKF,q_EKF,m_EKF,P_EKF,P_pose_prior,P_pose_posterior]=EKF_quick(N,delta_p,delta_q,y_mag,q_0,p_0,R_p,R_q,sigma_y,Lambda,Indices,sigma_lin,N_m,xl,xu,yl,yu,zl,zu);

p_RMSE_EKF=sqrt(mean(sum((p-p_EKF).^2,1)));
p_RMSE_DR=sqrt(mean(sum((p-p_DR).^2,1)));

mean_RMSEs_EKF(experiment,dataset)=p_RMSE_EKF;
mean_RMSEs_DR(experiment,dataset)=p_RMSE_DR;

end

disp(['Dataset ',num2str(dataset),', RMSE Algorithm 1: ',num2str(mean(mean_RMSEs_EKF(:,dataset)),'%4.2f'),'$\pm$',num2str(std(mean_RMSEs_EKF(:,dataset)),'%4.2f'),', RMSE dead reckoning: ',num2str(mean(mean_RMSEs_DR(:,dataset)),'%4.2f'),'$\pm$',num2str(std(mean_RMSEs_DR(:,dataset)),'%4.2f'),])

end
