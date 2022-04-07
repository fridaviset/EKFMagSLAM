function eta=quat2angleaxis(q)
%Copyright (C) 2022 by Frida Viset

coshalfalpha=q(1);
sinhalfalpha=norm(q(2:4));
halfalpha=atan2(sinhalfalpha,coshalfalpha);
alpha=halfalpha*2;
if sinhalfalpha<0.0001
    n=zeros(3,1);
else
    n=q(2:4)./sinhalfalpha;
end
if alpha>pi
    alpha=alpha-2*pi;
end
eta=alpha*n;
end