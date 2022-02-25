function [Indices, Lambda]=Lambda3D(N_m,uL1,lL1,uL2,lL2,uL3,lL3,sigma_SE,l_SE)
%Copyright (C) 2022 by Frida Viset

%Eigenvalues for the 1D squared exponential basis functions

max_test_val=25;
Indices=[];
for i=1:max_test_val
    for j=1:max_test_val
        for k=1:max_test_val
            Indices=[Indices; [i, j, k]];
        end
    end
end

j_1=Indices(:,1);
j_2=Indices(:,2);
j_3=Indices(:,3);

n=size(Indices,1);
eigv_candidates=zeros(n,1);

for j=1:n
    eigv=(pi*j_1(j)/(uL1-lL1)).^2+(pi*j_2(j)/(uL2-lL2)).^2+(pi*j_3(j)/(uL3-lL3)).^2;
    eigv_candidates(j)=spectral_density_SE(sqrt(eigv),sigma_SE,l_SE);
end

[eigvs_sorted, eigv_indices_sorted]=sort(eigv_candidates,'descend');
Indices=Indices(eigv_indices_sorted,:);
Indices=Indices(1:N_m,:);
    
Lambda=eye(N_m);
for j=1:N_m
    Lambda(j,j)=eigvs_sorted(j);
end

end