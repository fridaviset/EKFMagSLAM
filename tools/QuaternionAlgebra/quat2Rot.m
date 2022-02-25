function R=quat2Rot(q)
%Copyright (C) 2022 by Frida Viset


q0=q(1);
qv=q(2:4);
qvx=[0 -qv(3) qv(2);
     qv(3) 0 -qv(1);
     -qv(2) qv(1) 0];
R=qv*qv'+q0*q0*eye(3)+2*q0*qvx+qvx*qvx;


end