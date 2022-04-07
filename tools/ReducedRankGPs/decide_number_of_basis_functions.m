function N_m=decide_number_of_basis_functions(xl,xu,yl,yu,zl,zu,margin,sigma_SE,l_SE,sigma_y)
%Copyright (C) 2022 by Frida Viset

%Pick random locations to collect measurements and perform predictions
x_shift=xl+margin;
y_shift=yl+margin;
z_shift=zl+margin;
x_scale=(xu-xl)-2*margin;
y_scale=(yu-yl)-2*margin;
z_scale=(zu-zl)-2*margin;
N_meas=2000;
N_test=100;
x=diag([x_scale;y_scale;z_scale])*rand(3,N_meas+N_test)+[x_shift;y_shift;z_shift];
x_meas=x(:,1:N_meas);
x_test=x(:,N_meas+1:N_meas+N_test);

%Sample 100 test measurements in the first 100 locations
y_meas=mvnrnd(zeros(N_meas,1),SE_Kern(x_meas,x_meas,sigma_SE,l_SE)+sigma_y^2.*eye(N_meas));

%Use 100 test measurements to make prediction in remaining 100 locations
[mu, ~]=GaussianProcess(x_meas,y_meas,x_test,sigma_y,sigma_SE,l_SE);

%Start out with 10 basis functions, increase in increments of 10

for N_m=50:50:20000
    %Use the same 100 measurements to make a prediciton the same remaining 100 of the locations
    %with Reduced-rank approach
    [Indices, Lambda]=Lambda3D(N_m,xl,xu,yl,yu,zl,zu,sigma_SE,l_SE);
    %Reduced-rank GP
    phi_meas=Phi3D(x_meas,N_m,xl,xu,yl,yu,zl,zu,Indices);
    I=phi_meas'*phi_meas;
    iota=phi_meas'*y_meas';
    phi_test=Phi3D(x_test,N_m,xl,xu,yl,yu,zl,zu,Indices);
    mu_RR=phi_test*inv(I+sigma_y^2.*Lambda)*iota;
    
    %Evaluate the RMSE between the predictions
    RMSE=sqrt(mean((mu_RR-mu).^2));
    
    disp(['N_m: ',num2str(N_m),', RMSE:',num2str(RMSE)]);
    
    %If less than 10% of sigma_y, enough basis functions
    if RMSE<=sigma_y
        disp(['Sufficient amount of basis functions for current domain size and hyperparameters is estimated to be: ',num2str(N_m)]);
        break;
    end
end



end