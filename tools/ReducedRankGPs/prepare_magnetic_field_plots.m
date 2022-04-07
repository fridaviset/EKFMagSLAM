function [X,Y,Z,pointsVec,PhiVec,NablaPhi3D]=prepare_magnetic_field_plots(xl,xu,yl,yu,zl,zu,N_m,Indices,res,z)
%Copyright (C) 2022 by Frida Viset

x=xl:res:xu;
y=yl:res:yu;
[X,Y]=meshgrid(x,y);
Z=z*ones(size(X));
pointsVec=[X(:),Y(:),Z(:)]';
N_points=size(pointsVec,2);
PhiVec=[pointsVec', Phi3D(pointsVec,N_m,xl,xu,yl,yu,zl,zu,Indices)];
NablaPhi3D=cat(3,permute(repmat(eye(3),1,1,N_points),[3 1 2]), permute(Nabla_Phi3D(pointsVec,N_m,xl,xu,yl,yu,zl,zu,Indices),[3 1 2]));
end