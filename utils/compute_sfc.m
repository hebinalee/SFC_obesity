function sfc = compute_sfc(binconn, Nstep)
% 
% This function counts the number of all possible paths that connect different brain regions with specific step distances.
% The approach evaluates direct connections to indirect connections involving a varying number of step distances.
% 
% Inputs :  binconn - Binarized connectivity matrix (Nroi x Nroi)
% Outputs:  sfc     - Stepwise connectivity matrix (Nroi x Nroi x Nstep)
% 

[sfc, ~, ~] = findwalks(binconn);
sfc = sfc(:, :, 1:Nstep);
Nroi = size(binconn, 1);

for step = 1 : Nstep
    W = reshape(sfc(:,:,step), [Nroi, Nroi]);
    W = W - diag(diag(W));                           % off-diagonal
    sfc(:, :, step) = (W - mean(W(:))) / std(W(:));  % Normalization
end
end
