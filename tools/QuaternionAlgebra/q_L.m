function q_L_matrix=q_L(q)
%Copyright (C) 2022 by Frida Viset

qx=q(2);
qy=q(3);
qz=q(4);
OMEGA=[0 -qx -qy -qz;
    qx 0 -qz qy;
    qy qz 0 -qx;
    qz -qy qx 0];
q_L_matrix=(q(1)*eye(4)+OMEGA);
end