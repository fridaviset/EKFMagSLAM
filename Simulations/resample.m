function [new_q, new_p]=resample(old_q,old_p,weights)
%Copyright (C) 2022 by Frida Viset

M=length(weights);
indices=randsample(M,M,'true',weights);
new_q=old_q(:,:,indices);
new_p=old_p(:,:,indices);
end