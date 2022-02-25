function plot_projection_norm_magnetic_field(m,P_m,NablaPhi3D,X,Y,Z,fontsize)
%Copyright (C) 2022 by Frida Viset

N_points=size(NablaPhi3D,1);
N_m=size(m,1);
NablaPhi3DVec=reshape(NablaPhi3D,N_points*3,N_m);
FieldVec=NablaPhi3DVec*m;
FieldVec3=reshape(FieldVec,N_points,3);
N_points=size(FieldVec3,1);
Norm=zeros(N_points,1);
VarVec=zeros(N_points,1);
for i=1:N_points
    Norm(i)=norm(FieldVec3(i,:));
    NablaPhiHere=reshape(NablaPhi3D(i,:,:),3,N_m);
    VarVec(i)=trace(NablaPhiHere*P_m*NablaPhiHere');
end
Norm=reshape(Norm,size(X));
Var=reshape(VarVec,size(X));
surf(X,Y,Z,Norm,'AlphaData',1./Var,'FaceAlpha','flat','EdgeColor','none');
view(2);
xlabel('$x($m$)$','Interpreter','Latex','Fontsize',fontsize);
ylabel('$y($m$)$','Interpreter','Latex','Fontsize',fontsize);
zlabel('$z($m$)$','Interpreter','Latex','Fontsize',fontsize);
c=colorbar;
%title(c,'$\|H(p)\|_2^2$','Interpreter','Latex','Fontsize',fontsize);
end