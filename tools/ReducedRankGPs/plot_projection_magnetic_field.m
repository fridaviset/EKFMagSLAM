function plot_projection_magnetic_field(m,P_m,PhiVec,X,Y,Z)
%Copyright (C) 2022 by Frida Viset

FieldVec=PhiVec*m;
N=size(FieldVec,1);
VarVec=zeros(N,1);
for i=1:N
    VarVec(i)=PhiVec(i,:)*P_m*PhiVec(i,:)';
end
Field=reshape(FieldVec,size(X));
Var=reshape(VarVec,size(X));
surf(X,Y,Z,Field,'AlphaData',-Var,'FaceAlpha','flat','EdgeColor','none');
view(2);
colorbar;
end