function S_SE=spectral_density_SE(omega,sigma_SE,l_SE)
%Copyright (C) 2022 by Frida Viset

%This function is implemented based on
%Equation 21 in "Modeling and Interpolation of the Ambient Magnetic Field
%by Gaussian Processes" by Manon Kok and Arno Solin, published in
%IEEE TRANSACTIONS ON ROBOTICS, VOL. 34, NO. 4, AUGUST 2018
S_SE=sigma_SE^2*(2*pi*l_SE^2)^(3/2)*exp(-omega^2*l_SE^2/2);

end