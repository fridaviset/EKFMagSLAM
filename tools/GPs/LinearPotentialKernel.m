function K=LinearPotentialKernel(x,y,sigma_lin)
%Copyright (C) 2022 by Frida Viset

rows=size(x,2);
columns=size(y,2);
K=zeros(rows,columns);
for i=1:rows
    for j=1:columns
        K(i,j)=sigma_lin^2*(y(:,j)')*x(:,i);
    end
end
end