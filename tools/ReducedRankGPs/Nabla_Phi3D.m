function nabla_phi=Nabla_Phi3D(x,N_m,xl,xu,yl,yu,zl,zu,Indices)
%Copyright (C) 2022 by Frida Viset

%Basis functions for the SE kernel function
%Call this function first to find indices: [Indices, Lambda]=Lambda2D(m);
j1=Indices(:,1);
j2=Indices(:,2);
j3=Indices(:,3);

N=size(x,2);
nabla_phi=zeros(3,N_m,N);
for i=1:N
    nabla_phi(1,:,i)=1./sqrt(0.5*(xu-xl)).*cos(pi*j1'.*(x(1,i)-xl)./(xu-xl)).*pi.*j1'./(xu-xl)...
        ./sqrt(0.5*(yu-yl)).*sin(pi*j2'.*(x(2,i)-yl)./(yu-yl))...
        ./sqrt(0.5*(zu-zl)).*sin(pi*j3'.*(x(3,i)-zl)./(zu-zl));
    nabla_phi(2,:,i)=1./sqrt(0.5*(xu-xl)).*sin(pi*j1'.*(x(1,i)-xl)./(xu-xl))...
        ./sqrt(0.5*(yu-yl)).*cos(pi*j2'.*(x(2,i)-yl)./(yu-yl)).*pi.*j2'./(yu-yl)...
        ./sqrt(0.5*(zu-zl)).*sin(pi*j3'.*(x(3,i)-zl)./(zu-zl));
    nabla_phi(3,:,i)=1./sqrt(0.5*(xu-xl)).*sin(pi*j1'.*(x(1,i)-xl)./(xu-xl))...
        ./sqrt(0.5*(yu-yl)).*sin(pi*j2'.*(x(2,i)-yl)./(yu-yl))...
        ./sqrt(0.5*(zu-zl)).*cos(pi*j3'.*(x(3,i)-zl)./(zu-zl)).*pi.*j3'./(zu-zl);
end
end