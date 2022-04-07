%Copyright (C) 2022 by Frida Viset

close all; clear;

disp('Runnig algorithm 1 on first dataset from model ship with varying odometry noise');

load('Aligned_SquareLoop_1.mat');

seed=42;
rng(seed);

%Find room
margin=1;
xl=min(p(1,:))-margin;
xu=max(p(1,:))+margin;
yl=min(p(2,:))-margin;
yu=max(p(2,:))+margin;
zl=min(p(3,:))-margin;
zu=max(p(3,:))+margin;

%Number of basis functions used in Reduced-Rank approximation
N_m=50;

%Magnetic field params
sigma_SE=1;
l_SE=0.8;
sigma_lin=1;
sigma_y=0.1;

%Calculate Lambda and the order of indices used in the
%analytic basis functions of the Reduced-Rank Approximation
[Indices, Lambda]=Lambda3D(N_m,xl,xu,yl,yu,zl,zu,sigma_SE,l_SE);

%Experiment conditions
experiments=100;
sigma_xys=logspace(-4, -2, 10);
params=size(sigma_xys,2);

%Pre-allocate storage for MC results
p_Endpoint_EKF=zeros(params,experiments);
p_Endpoint_DR=zeros(params,experiments);
eta_Endpoint_DR=zeros(params,experiments);
eta_Endpoint_EKF=zeros(params,experiments);
norm_of_covar=zeros(params,experiments);

%Prepare folder for saving results in
time=clock;
foldername=['MCResults/Run',date,'-',num2str(time(4)),'-',num2str(time(5))];
mkdir(foldername);
mkdir([foldername,'/EPSs']);

for experiment=1:experiments
    for param=1:params
        
        %Set hyper-parameters
        R_p=0.0000001*eye(3);
        sigma_xy=sigma_xys(param);
        R_p(1:2,1:2)=sigma_xy*eye(2);
        R_q=0.00001*eye(3);
        bias=[0.003; 0.003; 0];
        
        %Simulate odometry
        [delta_p,delta_q,p_DR,q_DR,N]=pseudo_odometry(p,q,R_p,R_q,p_0,q_0,T,bias);
        
        %Adjust odometry noise to compensate for bias
        R_p(1:2,1:2)=(0.001+sigma_xy)*eye(2);
        
        %Prepare plotting
        fontsize=14;
        
        %Run the filter
        [p_EKF,q_EKF,m_EKF,P_EKF,P_pose_prior,P_pose_posterior]=EKF_quick(N,delta_p,delta_q,y_mag,q_0,p_0,R_p,R_q,sigma_y,Lambda,Indices,sigma_lin,N_m,xl,xu,yl,yu,zl,zu);
        
        %Find the max norm of the covariance
        P_pose_prior_norms=zeros(N,1);
        for t=1:N
            P_pose_prior_norms(t)=norm(P_pose_prior(:,:,t));
        end
        
        norm_of_covar(param,experiment)=max(P_pose_prior_norms);
        p_Endpoint_EKF(param,experiment)=norm(p(:,end)-p_EKF(:,end));
        p_Endpoint_DR(param,experiment)=norm(p(:,end)-p_DR(:,end));
        
        %Find the angle-axis error deviation
        err_DR=zeros(3,N);
        err_EKF=zeros(3,N);
        for t=1:N
            err_DR(:,t)=quat2angleaxis(quatprod(q_DR(:,t),[-q(1,t);q(2:4,t)]));
            err_EKF(:,t)=quat2angleaxis(quatprod(q_EKF(:,t),[-q(1,t);q(2:4,t)]));
        end
        
        eta_Endpoint_DR(param,experiment)=norm(err_DR(:,end));
        eta_Endpoint_EKF(param,experiment)=norm(err_EKF(:,end));
        
        disp(['experiment number: ',num2str(experiment),'/',num2str(experiments),', parameter number: ',num2str(param),'/',num2str(params)])
        
    end
end

save([foldername,'/Workspace.mat']);

figure; clf;
errorbar(sigma_xys,mean(p_Endpoint_DR,2),std(p_Endpoint_DR'),'b','linewidth',2);
hold on;
errorbar(sigma_xys,mean(p_Endpoint_EKF,2),std(p_Endpoint_EKF'),'k','linewidth',2);
xlabel('$\sigma_p^2$(m$^2$)','Interpreter','Latex','Fontsize',fontsize);
ylabel('$\|p_N-\hat{p}_N\|_2$(m)','Interpreter','Latex','Fontsize',fontsize);
legend({'Odometry','Algorithm 1'},'Interpreter','Latex','Fontsize',fontsize);
set(gca,'XScale','log');
ax = gca;
ax.FontSize = fontsize;
grid on;
saveas(gcf,[foldername,'/BoatPositionEstimationErrorsMCReps']);
saveas(gcf,[foldername,'/EPSs/BoatPositionEstimationErrorsMCReps'],'epsc');

figure; clf;
errorbar(sigma_xys,mean(eta_Endpoint_DR,2),std(eta_Endpoint_DR'),'b','linewidth',2);
hold on;
errorbar(sigma_xys,mean(eta_Endpoint_EKF,2),std(eta_Endpoint_EKF'),'k','linewidth',2);
ax = gca;
ax.FontSize = fontsize;
grid on;
xlabel('$\sigma_p^2$(m$^2$)','Interpreter','Latex','Fontsize',fontsize);
ylabel('$\eta_N$','Interpreter','Latex','Fontsize',fontsize);
legend({'Odometry','Algorithm 1'},'Interpreter','Latex','Fontsize',fontsize);
set(gca,'XScale','log');
saveas(gcf,[foldername,'/BoatOrientationEstimationErrors']);
saveas(gcf,[foldername,'/EPSs/BoatOrientationEstimationErrorsMCReps'],'epsc');

figure; clf;
errorbar(sigma_xys,mean(norm_of_covar,2),std(norm_of_covar'),'k','linewidth',2);
ax = gca;
ax.FontSize = fontsize;
grid on;
xlabel('$\sigma_p^2$(m$^2$)','Interpreter','Latex','Fontsize',fontsize);
ylabel('Max norm of preditictive covariance','Interpreter','Latex','Fontsize',fontsize);
set(gca,'XScale','log');
saveas(gcf,[foldername,'/BoatPriorCovarianceNorm']);
saveas(gcf,[foldername,'/EPSs/BoatPriorCovarianceNorm'],'epsc');
