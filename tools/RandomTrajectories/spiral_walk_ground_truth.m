function [p,p_0,q,q_0,N]=spiral_walk_ground_truth(T,xl,xu,yl,yu,margin,z,v,gap)
%Copyright (C) 2022 by Frida Viset

%xl - lower limit x, a.k.a. the position of the left wall [m]
%xu - upper limit x, a.k.a. the position of the right wall [m]
%yl - lower limit y, a.k.a. the position of the lower wall [m]
%yu - uppper limit y, a.k.a. the position of the upper wall [m]
%margin - distance from trajectory to wall [m]
%v - walking speed [m/s]
%z - position along z-axis of trajectory [m]
%gap - distance to the next loop of the spiral. Needs to be larger than
%v*T for the spiral to make sense.

%length - Length of trajectory in meter
length=xu-xl-2*margin;

%height - %Height of trajectory in meter
height=yu-yl-2*margin;

%Distance traversed per sample
diff_p=v*T;

%Make the four edges of the rectangle
p1=[diff_p:diff_p:(length); zeros(1,floor(length/diff_p)); zeros(1,floor(length/diff_p))];
p2=[zeros(1,floor(height/diff_p)); diff_p:diff_p:(height); zeros(1,floor(height/diff_p))]+p1(:,end);
p3=[-(diff_p:diff_p:(length)); zeros(1,floor(length/diff_p)); zeros(1,floor(length/diff_p))]+p2(:,end);

%Concatenate the edges to make first part of spiral
p=[p1 p2 p3];

for i=1:floor(min((length)./(2*gap),(height)./(2*gap)))
    length=length-gap;
    height=height-gap;
    p1=[zeros(1,floor(height/diff_p)); -(diff_p:diff_p:(height)); zeros(1,floor(height/diff_p))]+p(:,end);
    p2=[diff_p:diff_p:(length); zeros(1,floor(length/diff_p)); zeros(1,floor(length/diff_p))]+p1(:,end);
    length=length-gap;
    height=height-gap;
    p3=[zeros(1,floor(height/diff_p)); diff_p:diff_p:(height); zeros(1,floor(height/diff_p))]+p2(:,end);
    p4=[-(diff_p:diff_p:(length)); zeros(1,floor(length/diff_p)); zeros(1,floor(length/diff_p))]+p3(:,end);
    p=[p p1 p2 p3 p4];
end

%Shift all positions according to the initial position
p=p+[xl+margin; yl+margin; z];
p_0=p(:,1);

%Find the length of the trajectory
N=size(p,2);

%Create very boring orientation odometry (constant orientation)
q_0=[1; 0; 0; 0];
q(:,1)=q_0;
for t=2:N
    q(:,t)=q(:,t-1);
end
end