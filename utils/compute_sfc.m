function sfc = compute_sfc(binconn, Nstep)
[sfc, ~, ~] = findwalks(binconn);
sfc = sfc(:, :, 1:Nstep);
Nroi = size(binconn, 1);

for step = 1 : Nstep
    W = reshape(sfc(:,:,step), [Nroi, Nroi]);
    W = W - diag(diag(W));                           % off-diagonal
    sfc(:, :, step) = (W - mean(W(:))) / std(W(:));  % Normalization
end
end
