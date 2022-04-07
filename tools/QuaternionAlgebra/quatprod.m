function prod=quatprod(p,q)
%Copyright (C) 2022 by Frida Viset

q0=q(1);
p0=p(1);
qv=q(2:4);
pv=p(2:4);
pvx=[0 -pv(3) pv(2);
     pv(3) 0 -pv(1);
     -pv(2) pv(1) 0];
prod=[p0*q0-pv'*qv;
      p0*qv+q0*pv+pvx*qv];
end