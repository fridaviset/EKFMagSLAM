function [y_gyr,y_acc,y_mag]=simulate_IMU(N,q,acc,omega,m,g,R_gyr,R_acc,R_mag)
%Copyright (C) 2022 by Frida Viset

%%Inputs
%N - number of time steps
%q - orientation, 4 times N
%acc - acceleration, 3 times N
%omega - angular velocity, 3 times N
%m - magnetic field, 3 times 1
%g - gravity field, 3 times 1
%R_gyr - variance of gyroscope measurements, 3 times 3
%R_acc - variance of accelerometer measurement, 3 times 3
%R_mag - variance of magnetometer measurements, 3 times 3

%%Outputs
%y_gyr - simulated gyroscope measurements, 3 times N
%y_acc - simulated accelerometer measurements, 3 times N
%y_mag - simulated magnetometer measurements, 3 times N

%Assigne meas vectors
y_acc=zeros(3,N);
y_gyr=zeros(3,N);
y_mag=zeros(3,N);

%Simulate measurements
for t=1:N
y_acc(:,t)=quat2Rot(q(:,t))'*(acc(:,t)-g)+mvnrnd(zeros(3,1),R_acc)';
y_gyr(:,t)=quat2Rot(q(:,t))'*omega(:,t)+mvnrnd(zeros(3,1),R_gyr)';
y_mag(:,t)=quat2Rot(q(:,t))'*(m)+mvnrnd(zeros(3,1),R_mag)';
end


end