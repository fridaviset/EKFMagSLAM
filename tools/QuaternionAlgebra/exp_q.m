function q=exp_q(eta)
%Copyright (C) 2022 by Frida Viset

q=zeros(4,1);
alpha=norm(eta);
if alpha<0.0001
q(1)=1;
else
q(1)=cos(norm(eta)./2);
q(2:4)=sin(norm(eta)./2)*eta./norm(eta);
end
end