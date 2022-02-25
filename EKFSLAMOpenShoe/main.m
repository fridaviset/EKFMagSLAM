clear;
close all;

%Copyright (C) 2022 by Frida Viset

disp('Running Algorithm 1 on data from foot-mounted sensor');

time=clock;

%Set room dimensions
margin=10;

%Load dataset
load('OpenShoeOdometryAndMagField.mat');

%Run dead reckoning on the simplified odometry
[p_DR,q_DR]=Simple_Dead_Reckoning(delta_p,delta_q,p_0,q_0);

%Estimate the dimensions of the problem based on the first round
traj=p_DR(:,1:800);
xl=min(traj(1,:))-margin;
xu=max(traj(1,:))+margin;
yl=min(traj(2,:))-margin;
yu=max(traj(2,:))+margin;
zl=min(traj(3,:))-margin;
zu=max(traj(3,:))+margin;

%Number of basis functions used in Reduced-Rank approximation
N_m=2000;

%Magnetic field params
sigma_SE=1;
l_SE=2;
sigma_lin=1;
sigma_y=0.1;

%Calculate Lambda and the order of indices used in the
%analytic basis functions of the Reduced-Rank Approximation
[Indices, Lambda]=Lambda3D(N_m,xl,xu,yl,yu,zl,zu,sigma_SE,l_SE);

%Noise parameters
R_p=0.001*eye(3)*T;
R_q=0.00001*eye(3)*T;

%Init covariances
P_0=0.0001*eye(6);

%Prepare folder for saving results in
foldername=['Results/Run',date,'-',num2str(time(4)),'-',num2str(time(5))];
mkdir(foldername);
mkdir([foldername,'/EPSs']);
mkdir([foldername,'/JPGs']);
mkdir([foldername,'/Snapshots']);

%Prepare plotting
fontsize=14;

%Run the filter
[p_EKF,q_EKF,m,P,P_pose_prior,P_pose_posterior]=EKF_Fast(N,T,foldername,delta_p,delta_q,y_mag,q_0,p_0,R_p,R_q,sigma_y,Lambda,Indices,sigma_lin,N_m,xl,xu,yl,yu,zl,zu,fontsize);

close all;

figure; clf;
hold on;
plot3(p_EKF(1,:),p_EKF(2,:),p_EKF(3,:),'k');
plot3(p_DR(1,:),p_DR(2,:),p_DR(3,:),'b');
legend({'Algorithm 1','Odometry'},'Interpreter','Latex','Fontsize',fontsize)
view(2);
xlabel('$x($m$)$','Interpreter','Latex','Fontsize',fontsize);
ylabel('$y($m$)$','Interpreter','Latex','Fontsize',fontsize);
zlabel('$z$(m$)$','Interpreter','Latex','Fontsize',fontsize);
axis equal;
xlim([xl xu]);
ylim([yl yu]);
saveas(gcf,[foldername,'/FootTrajectoryComparisons']);
filename=[foldername,'/EPSs/FootTrajectoryComparisons.eps'];
exportgraphics(gca,filename,'BackgroundColor','none','ContentType','vector');
filename=[foldername,'/JPGs/FootTrajectoryComparisons.jpg'];
exportgraphics(gca,filename,'BackgroundColor','none','ContentType','vector');

%Pre-generate some plotting matrices
res=0.25; zs=linspace(min(p_EKF(3,:)),max(p_EKF(3,:)),10);
P_m=P;
figure; clf;
for i=1:length(zs)
[X,Y,Z,pointsVec,PhiVec,NablaPhi3D]=prepare_magnetic_field_plots(xl,xu,yl,yu,zl,zu,N_m,Indices,res,zs(i));
plot_projection_norm_magnetic_field(m,P_m,NablaPhi3D,X,Y,Z,fontsize);
hold on;
end
plot3(p_EKF(1,:),p_EKF(2,:),p_EKF(3,:),'k');
plot3(p_EKF(1,:),p_EKF(2,:),0.*p_EKF(3,:)+3,'k');
view(2);
axis equal;
xlim([xl xu]);
ylim([yl yu]);
cmin=0;
cmax=1;
caxis([cmin cmax]);
grid off;
saveas(gcf,[foldername,'/MagFieldEstFoot']);
filename=[foldername,'/EPSs/MagFieldEstFoot.eps'];
exportgraphics(gca,filename,'BackgroundColor','none','ContentType','vector');
filename=[foldername,'/EPSs/MagFieldEstFoot.pdf'];
exportgraphics(gca,filename,'BackgroundColor','none','ContentType','vector');
filename=[foldername,'/JPGs/MagFieldEstFoot.jpg'];
exportgraphics(gca,filename,'BackgroundColor','none','ContentType','vector');
