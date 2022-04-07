function K=CurlFreeKernel(x,y,sigma_SE,l_SE)
%Copyright (C) 2022 by Frida Viset

rows=size(x,2);
columns=size(y,2);
K=zeros(3*rows,3*columns);
for i=1:rows
    for j=1:columns
        p=x(:,i)-y(:,j);
        nom=norm(p)^2;
        den=2*l_SE^2;
        core=nom/den;
        K(3*(i-1)+1:3*(i-1)+3,3*(j-1)+1:3*(j-1)+3)=sigma_SE^2./l_SE^2.*exp(-core)*(eye(3)-p*p'./l_SE^2);
    end
end
end