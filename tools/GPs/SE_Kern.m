function K=SE_Kern(x,y,sigma_SE,l_SE)
%Copyright (C) 2022 by Frida Viset

rows=size(x,2);
columns=size(y,2);
K=zeros(rows,columns);
for i=1:rows
    for j=1:columns
        p=x(:,i)-y(:,j);
        nom=norm(p)^2;
        den=2*l_SE^2;
        core=nom/den;
        K(i,j)=sigma_SE^2*exp(-core);
    end
end
end