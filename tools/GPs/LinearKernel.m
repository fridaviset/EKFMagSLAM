function K=LinearKernel(x,y,sigma_lin)
%Copyright (C) 2022 by Frida Viset

rows=size(x,2);
columns=size(y,2);
K=zeros(3*rows,3*columns);
for i=1:rows
    for j=1:columns
        K(3*(i-1)+1:3*(i-1)+3,3*(j-1)+1:3*(j-1)+3)=sigma_lin^2*eye(3);
    end
end
end