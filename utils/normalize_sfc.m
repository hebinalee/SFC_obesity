function norm_sfc = normalize_sfc(sfc)
norm_sfc = sfc;
steps = size(norm_sfc,3);
for i = 1 : steps
	stepfc = norm_sfc(:,:,i);
    norm_sfc(:,:,i) = (stepfc - mean(stepfc(:))) / std(stepfc(:));
end