%Copyright (C) 2022 by Frida Viset
clear;
close all;

disp('Running simulations');
tic;

seed=42;
rng(seed);

%Prepare folder for saving results in
time=clock;
foldername=['LocalisationMCResults/Run',date,'-',num2str(time(4)),'-',num2str(time(5))];
mkdir(foldername);
mkdir([foldername,'/EPSs']);
mkdir([foldername,'/InitialDistributions']);
mkdir([foldername,'/InitialDistributions/EPSs']);
mkdir([foldername,'/Trajectories'])
mkdir([foldername,'/PositionEstimationErrors'])
mkdir([foldername,'/OrientationEstimationErrors'])

%lengthscales=logspace(log(0.03)./log(10),log(1.5)./log(10),10);
%params=size(lengthscales,2);
initial_errors=[0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5];
params=size(initial_errors,2);
N_p_s=[100 200 500];
particle_indices=size(N_p_s,2);
experiments=100;

p_Endpoint_EKF=zeros(params,experiments);
p_Endpoint_DR=zeros(params,experiments);
p_Endpoint_PF=zeros(params,experiments);

%Set room dimensions
margin=0.5;
xl=-1.5;
xu=1.5;
yl=-1.5;
yu=1.5;
zl=-0.5;
zu=0.5;
T=0.1;

%Simulate trajectory instead
z=0.0; v=0.5; gap=0.1;
[p,p_0_gt,q,q_0,N]=rectangular_walk_ground_truth(T,xl,xu,yl,yu,margin,z,v);

%Decide upon the number of basis functions that causes the GP
%approximation error of 100 samples to be less than the measurement noise
%Magnetic field params
sigma_SE=0.1;
sigma_y=0.03;
sigma_lin=0;
l_SE=0.3;
N_m=decide_number_of_basis_functions(xl,xu,yl,yu,zl,zu,margin,sigma_SE,l_SE,sigma_y);
[Indices, Lambda]=Lambda3D(N_m,xl,xu,yl,yu,zl,zu,sigma_SE,l_SE);

for param=1:params
    
    close all;
    initial_error=initial_errors(param);
    mkdir([foldername,'/InitialDistributions/param',num2str(param)]);
    mkdir([foldername,'/InitialDistributions/EPSs/param',num2str(param)]);
    mkdir([foldername,'/Trajectories/param',num2str(param)]);
    mkdir([foldername,'/PositionEstimationErrors/param',num2str(param)]);
    mkdir([foldername,'/OrientationEstimationErrors/param',num2str(param)]);
    
    for experiment=1:experiments
        
        %Simulate magnetic field
        P_0=[sigma_lin^2*eye(3), zeros(3,N_m);
            zeros(N_m,3), Lambda];
        m_sim=mvnrnd(zeros(N_m+3,1),P_0)';
        
        %Simulate magnetic field measurements
        y_mag=zeros(3,N);
        for t=1:N
            NablaPhi=[eye(3),Nabla_Phi3D(p(:,t),N_m,xl,xu,yl,yu,zl,zu,Indices)];
            f=NablaPhi*m_sim;
            y_mag(:,t)=quat2Rot(q(:,t))'*f+mvnrnd(zeros(3,1),sigma_y^2);
        end
        
        %Number of particles used in visuals
        M=N_p_s(end);
        
        %Noise parameters
        R_p=0.001*eye(3)*T;
        R_q=0.0001*eye(3)*T;
        
        %Init covariances
        P_0=0.001*eye(6);
        P_0(1:2,1:2)=initial_error*eye(2);
        
        %Add error to the position estimate
        p_0=p_0_gt+[sqrt(0.5)*initial_error; sqrt(0.5)*initial_error; 0];
        
        %Find odometry measurements
        bias=[0; 0; 0];
        [delta_p,delta_q,p_DR,q_DR,N]=pseudo_odometry(p,q,R_p,R_q,p_0,q_0,T,bias);
        
        %Final plots for result presentation
        fontsize=14;
        
        %Run filter with higher measurement error to better visualise the
        %influence of the measurement on the PF distribution
        sigma_y=0.3;
        
        %Run the filters
        [p_PF,q_PF,p_EKF,q_EKF]=PF_and_EKF_localisation(N,T,foldername,delta_p,delta_q,y_mag,q_0,p_0,P_0,p,R_p,R_q,m_sim,sigma_y,Indices,N_m,xl,xu,yl,yu,zl,zu,M,fontsize,experiment,param);
        
        %Repeat for all number of particles
        for particle_index=1:particle_indices
            M=N_p_s(particle_index);
            [p_PF,q_PF,p_EKF,q_EKF]=PF_and_EKF_localisation(N,T,foldername,delta_p,delta_q,y_mag,q_0,p_0,P_0,p,R_p,R_q,m_sim,sigma_y,Indices,N_m,xl,xu,yl,yu,zl,zu,M,fontsize,2,param);
            p_RMSE_PF=sqrt(mean(sum((p-p_PF).^2,1)));
            p_Endpoint_PF(param,experiment,particle_index)=p_RMSE_PF(end);
        end
        
        close all;
        
        if experiment==1
            
            figure; clf;
            title(['Predictive position error $',num2str(initial_error),'$ meters'],'Interpreter','Latex');
            hold on;
            plot3(p_EKF(1,:),p_EKF(2,:),p_EKF(3,:),'k');
            plot3(p_PF(1,:),p_PF(2,:),p_PF(3,:),'g');
            plot3(p_DR(1,:),p_DR(2,:),p_DR(3,:),'b');
            plot3(p(1,:),p(2,:),p(3,:),'r');
            view(2);
            xlabel('x/[m]','Interpreter','Latex');
            ylabel('y/[m]','Interpreter','Latex');
            zlabel('z/[m]','Interpreter','Latex');
            legend({'EKF','PF','Dead Reckoning'},'Interpreter','Latex','Fontsize',fontsize);
            saveas(gcf,[foldername,'/Trajectories/param',num2str(param),'/experiment',num2str(experiment)]);
            
            figure; clf;
            plot(sqrt(sum((p-p_EKF).^2,1)),'k');
            hold on;
            plot(sqrt(sum((p-p_PF).^2,1)),'g');
            plot(sqrt(sum((p-p_DR).^2,1)),'b');
            xlabel('$t$(s)','Interpreter','Latex','Fontsize',fontsize);
            ylabel('$\|p_t-\hat{p}_t\|_2$','Interpreter','Latex','Fontsize',fontsize);
            legend({'EKF','PF','Dead Reckoning'},'Interpreter','Latex','Fontsize',fontsize);
            saveas(gcf,[foldername,'/PositionEstimationErrors/param',num2str(param),'/experiment',num2str(experiment)]);
            
            %Find the angle-axis error deviation
            err_DR=zeros(3,N);
            err_EKF=zeros(3,N);
            err_PF=zeros(3,N);
            for t=1:N
                err_DR(:,t)=quat2angleaxis(quatprod(q_DR(:,t),[-q(1,t);q(2:4,t)]));
                err_EKF(:,t)=quat2angleaxis(quatprod(q_EKF(:,t),[-q(1,t);q(2:4,t)]));
                err_PF(:,t)=quat2angleaxis(quatprod(q_PF(:,t),[-q(1,t);q(2:4,t)]));
            end
            
            figure; clf;
            subplot(3,1,1);
            plot(err_DR');
            legend({'roll','pitch','yaw'},'Interpreter','Latex','Fontsize',fontsize);
            title('Dead Reckoning','Interpreter','Latex','Fontsize',fontsize);
            xlabel('$t$(s)','Interpreter','Latex','Fontsize',fontsize);
            ylabel('$\eta_t$','Interpreter','Latex','Fontsize',fontsize);
            subplot(3,1,2);
            plot(err_PF');
            legend({'roll','pitch','yaw'},'Interpreter','Latex','Fontsize',fontsize);
            title('PF','Interpreter','Latex','Fontsize',fontsize);
            xlabel('$t$(s)','Interpreter','Latex','Fontsize',fontsize);
            ylabel('$\eta_t$','Interpreter','Latex','Fontsize',fontsize);
            subplot(3,1,3);
            plot(err_EKF');
            legend({'roll','pitch','yaw'},'Interpreter','Latex','Fontsize',fontsize);
            title('EKF','Interpreter','Latex','Fontsize',fontsize);
            xlabel('$t$(s)','Interpreter','Latex','Fontsize',fontsize);
            ylabel('$\eta_t$','Interpreter','Latex','Fontsize',fontsize);
            saveas(gcf,[foldername,'/OrientationEstimationErrors/param',num2str(param),'/experiment',num2str(experiment)]);
            
        end
        
        p_RMSE_EKF=sqrt(mean(sum((p-p_EKF).^2,1)));
        p_RMSE_DR=sqrt(mean(sum((p-p_DR).^2,1)));
        
        p_Endpoint_EKF(param,experiment)=p_RMSE_EKF(end);
        p_Endpoint_DR(param,experiment)=p_RMSE_DR(end);
        
        
        eta_Endpoint_EKF(param,experiment)=norm(err_EKF(:,end));
        eta_Endpoint_DR(param,experiment)=norm(err_DR(:,end));
        eta_Endpoint_PF(param,experiment)=norm(err_PF(:,end));
        
        disp(['experiment number: ',num2str(experiment),'/',num2str(experiments),', parameter number: ',num2str(param),'/',num2str(params)])
        
    end
    
end

figure; clf;
errorbar(initial_errors,mean(p_Endpoint_DR,2),std(p_Endpoint_DR'),'b','linewidth',2);
hold on;
errorbar(initial_errors,mean(p_Endpoint_EKF,2),std(p_Endpoint_EKF'),'k','linewidth',2);
color1=[220, 38, 127]./255;
color2=[255, 176, 0]./255;
color3=[100, 143, 255]./255;
colors=[color1; color2; color3];
labels={'Odometry','Algorithm 1'};
for particle_index=1:particle_indices
    labels{particle_index+2}=['PF, $N_p=$',num2str(N_p_s(particle_index))];
    color=colors(particle_index,:);
    errorbar(initial_errors,mean(p_Endpoint_PF(:,:,particle_index),2),std(p_Endpoint_PF(:,:,particle_index)'),'Color',color,'linewidth',2);
end
xlabel('Predictive position estimation error (m)','Interpreter','Latex','Fontsize',fontsize);
ylabel('$\|p_N-\hat{p}_N\|_2$(m)','Interpreter','Latex','Fontsize',fontsize);
legend(labels,'Interpreter','Latex','Fontsize',fontsize,'Location','northwest');
ax = gca;
ax.FontSize = fontsize;
grid on;
saveas(gcf,[foldername,'/LocPositionEstimationErrorsMCReps']);
saveas(gcf,[foldername,'/EPSs/LocPositionEstimationErrorsMCReps'],'epsc');
saveas(gcf,[foldername,'/LocPositionEstimationErrorsMCReps'],'png');

timespent=toc;
save('timespent.mat','timespent');
save([foldername,'/Workspace.mat']);

