function [quat, g]=my_init_orientation(u,calib_time)
%Copyright (C) 2022 by Frida Viset

%This function is an implementation of Theorem 4.1 in Jeroen Hoel's PhD
%thesis "Sensor Fusion and Calibration of Inertial Sensors, Vision,
%Ultra-Wideband and GPS"

A=zeros(4);
%Assume the acceleration vector is only caused by gravity
for t=1:calib_time
    g_w=[0; 0; 0; norm(u(1:3,t))];
    g_b=[0; u(1:3,t)];
    A=A+q_L(g_b)*q_R(g_w);
end
[V,~] = eig(A);
quat=V(:,1);

%Find an estimate of the gravity magnitude
g=[0; 0; norm(u(1:3,t))];
end