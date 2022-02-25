function y_mag=simulate_y_mag(p,q,sigma_SE,sigma_lin,l_SE,sigma_y)
%Copyright (C) 2022 by Frida Viset

N=size(p,2);

%Simulate magnetic field measurements
K_SE=CurlFreeKernel(p,p,sigma_SE,l_SE);
K_l=LinearKernel(p,p,sigma_lin);
K_noise=sigma_y^2*eye(N*3);
K=K_SE+K_l+K_noise;
y_mag_vectorized=mvnrnd(zeros(N*3,1),K);
y_mag=reshape(y_mag_vectorized,[3,N]);

%And move the simulated magnetic field to body frame
for t=1:N
    y_mag(:,t)=quat2Rot(q(:,t))'*y_mag(:,t);
end

end