%Copyright (C) 2022 by Frida Viset

close all; clear;

load('Aligned_SquareLoop_1.mat');

disp('Runnig algorithm 1 on first dataset from model ship with detailed statistics');

seed=42;
rng(seed);

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

%Prepare folder for saving results in
time=clock;
foldername=['Results/Run',date,'-',num2str(time(4)),'-',num2str(time(5))];
mkdir(foldername);
mkdir([foldername,'/Snapshots']);
mkdir([foldername,'/EPSs']);
mkdir([foldername,'/JPGs']);

%Prepare plotting
fontsize=14;

%Run the filter
[p_EKF,q_EKF,m_EKF,P_EKF,P_pose_prior,P_pose_posterior]=EKF_quick(N,delta_p,delta_q,y_mag,q_0,p_0,R_p,R_q,sigma_y,Lambda,Indices,sigma_lin,N_m,xl,xu,yl,yu,zl,zu);

close all;

figure; clf;
hold on;
plot3(p_EKF(1,:),p_EKF(2,:),p_EKF(3,:),'k');
plot3(p_DR(1,:),p_DR(2,:),p_DR(3,:),'b');
plot3(p(1,:),p(2,:),p(3,:),'r');
legend({'Algorithm 1','Dead reckoning','Ground truth'},'Interpreter','Latex','Fontsize',fontsize);
view(2);
xlabel('$x($m$)$','Interpreter','Latex','Fontsize',fontsize);
ylabel('$y($m$)$','Interpreter','Latex','Fontsize',fontsize);
zlabel('$z($m$)$','Interpreter','Latex');
axis equal;
saveas(gcf,[foldername,'/TrajectoryComparisons']);
saveas(gcf,[foldername,'/EPSs/TrajectoryComparisons'],'epsc');
saveas(gcf,[foldername,'/JPGs/TrajectoryComparisons'],'jpg');

figure; clf;
time=T:T:T*N;
plot(time,sqrt(sum((p-p_EKF).^2,1)),'k');
hold on;
plot(time,sqrt(sum((p-p_DR).^2,1)),'b');
xlabel('time(s)','Interpreter','Latex','Fontsize',fontsize);
ylabel('$\|p_t-\hat{p}_t\|_2$(m)','Interpreter','Latex','Fontsize',fontsize);
legend({'Algorithm 1','Dead reckoning'},'Interpreter','Latex','Fontsize',fontsize);
saveas(gcf,[foldername,'/PositionEstimationErrors']);
saveas(gcf,[foldername,'/EPSs/PositionEstimationErrors'],'epsc');

p_RMSE_EKF=sqrt(mean(sum((p-p_EKF).^2,1)));
p_RMSE_DR=sqrt(mean(sum((p-p_DR).^2,1)));

disp(['RMSE EKF: ',num2str(p_RMSE_EKF),', RMSE DR: ',num2str(p_RMSE_DR)]);

%Pre-generate some plotting matrices
res=0.01; z=-0.2;
P_m=P_EKF(7:end,7:end);
[X,Y,Z,pointsVec,PhiVec,NablaPhi3D]=prepare_magnetic_field_plots(xl,xu,yl,yu,zl,zu,N_m,Indices,res,z);
figure; clf;
plot_projection_norm_magnetic_field(m_EKF,P_m,NablaPhi3D,X,Y,Z,fontsize);
hold on;
plot3(p_EKF(1,:),p_EKF(2,:),p_EKF(3,:)+50,'k');
axis equal;
view(2);
caxis([0.4 0.9]);
xlim([xl xu]);
ylim([yl yu]);
filename=[foldername,'/EPSs/MagFieldEstBoat.eps'];
exportgraphics(gca,filename,'BackgroundColor','none','ContentType','vector');
filename=[foldername,'/EPSs/MagFieldEstBoat.jpg'];
exportgraphics(gca,filename,'BackgroundColor','none','ContentType','vector');

%Compare with the real norm
figure; clf; 
scatter3(p(1,:),p(2,:),p(3,:),10+0.*y_mag_norm,y_mag_norm,'filled');
xlabel('$x($m$)$','Interpreter','Latex');
ylabel('$y($m$)$','Interpreter','Latex');
zlabel('$z($m$)$','Interpreter','Latex');
c=colorbar;
caxis([0.4 0.9]);
axis equal; 
xlim([xl xu]);
ylim([yl yu]);
view(2);
filename=[foldername,'/EPSs/MagFieldMeasBoat.eps'];
exportgraphics(gca,filename,'BackgroundColor','none','ContentType','vector');