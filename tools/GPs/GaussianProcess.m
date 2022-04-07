function [mu, var]=GaussianProcess(x,y,x_s,sigma_y,sigma_SE,l_SE)
%Copyright (C) 2022 by Frida Viset

N=size(x,2);
K=SE_Kern(x,x,sigma_SE,l_SE);
Ks=SE_Kern(x_s,x,sigma_SE,l_SE);
Kss=SE_Kern(x_s,x_s,sigma_SE,l_SE);
mu=Ks*inv(K+sigma_y^2.*eye(N))*y';
var=Kss-Ks*inv(K+sigma_y^2.*eye(N))*Ks';

end


