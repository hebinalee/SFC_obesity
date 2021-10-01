function high_low_risk()
inpath = 'X:/path/myfolder/inputs/';
outpath = 'X:/path/myfolder/outputs/';
Nsub = 301;
Nroi = 246;
Nstep = 7;
Nperm = 1000;


%% 1) Divide participants into two group based on WHR
disp(['## Divide subjects into low/high-health risk groups - processing', newline]);
load([inpath, 'a_dataset.mat'])
group = (sex==1 & obesity(:,2)<0.95) | (sex==2 & obesity(:,2)<0.8);             % low-risk
group = group + 2*((sex==1 & obesity(:,2)>1.0) | (sex==2 & obesity(:,2)>0.86)); % high-risk
group1_idx = find(group == 1);
group2_idx = find(group == 2);
Nsubsample = length(group2_idx);
save([outpath, 'high_low/group.mat'], 'group')


%% 2) Subject bootstrapping
disp(['## Subject bootstrapping for low-risk group - processing', newline]);
subsample = zeros(Nsubsample, Nperm);
bootstrap_group = zeros(Nsub, Nperm);
for i = 1 : Nperm
    subsample(:, i) = randperm(length(group1_idx), Nsubsample);
    if length(unique(subsample(:, i))) ~= Nsubsample
        disp(['i = ', num2str(i), ' : Repetition error during sampling'])
    end
    bootstrap_group(group1_idx(subsample(:,i)), i) = 1;
    bootstrap_group(group2_idx, i) = 2;
end
save([outpath, 'high_low/bootstrap_subsamples.mat'], 'subsample', 'bootstrap_group')


%% 3) Compute group average SFC matrix
disp(['## Compute group average SFC matrix during iterations - processing', newline]);
for i = 1 : Nperm
    if rem(i,100) == 0
        disp(['i = ', num2str(i)]);
    end
    seed_idx = find(bootstrap_selected(:,i)==1);    
    grpmean_SFC = cell(2,1);
    grpmean_SFC{1} = zeros(length(seed_idx), Nroi, Nstep);
    grpmean_SFC{2} = zeros(length(seed_idx), Nroi, Nstep);
    for j = 1 : Nsubsample
        if group(subsample(j,i)) == 1
            load([inpath, '2.sfc5/sub', pad(num2str(subsample(j,i), '%d'), 3, 'left', '0'), '.mat']);
            grpmean_SFC{1} = grpmean_SFC{1} + normalize_sfc(sfc(seed_idx,:,1:Nstep));
        elseif group(subsample(j,i)) == 2
            load([inpath, '2.sfc5/sub', pad(num2str(subsample(j,i), '%d'), 3, 'left', '0'), '.mat']);
            grpmean_SFC{2} = grpmean_SFC{2} + normalize_sfc(sfc(seed_idx,:,1:Nstep));
        end
    end
    grpmean_SFC{1} = grpmean_SFC{1} / sum(group(subsample(:,i)) == 1);
    grpmean_SFC{2} = grpmean_SFC{2} / sum(group(subsample(:,i)) == 2);
    save([outpath, 'high_low/sfc/iteration_', num2str(i), '.mat'], 'grpmean_SFC')
end


%% 4) Randomly pick representative iteration
disp(['## Randomly pick representative iteration - processing', newline]);
i = randi(Nperm, 1, 1);
load([inpath, 'a_dataset.mat'])
load([outpath, 'high_low/sfc/iteration_', num2str(i), '.mat'])
grpmean_DC(1,:,:) = mean(grpmean_SFC{1}, 1);
grpmean_DC(2,:,:) = mean(grpmean_SFC{2}, 1);
group = bootstrap_group(:, i);


%% 5) Group difference test: ROI-level
disp(['## Group difference test: ROI-level - processing', newline]);
load([outpath, 'wholesub_ROI_dc.mat']);
Nroi = 224;
H = zeros(Nroi, Nstep);
P = zeros(Nroi, Nstep);
T = zeros(Nroi, Nstep);

for step = 1 : Nstep
    ob_dc = [roi_dc(group==2,1:210,step), mean_subcortical(roi_dc(group==2,:,step)')'];
    hw_dc = [roi_dc(group==1,1:210,step), mean_subcortical(roi_dc(group==1,:,step)')'];
    for roi = 1 : Nroi
        [~,p,~,stats] = ttest2(ob_dc(:,roi), hw_dc(:,roi));
        P(roi, step) = p;
        T(roi, step) = stats.tstat;
    end
    [selected, ~, ~, corrected_P] = fdr_bh(P(:, step), 0.05);
    H(:, step) = selected;
    P(:, step) = corrected_P;
end
save([outpath, 'high_low/groupdiff_ROI_ttest.mat'], 'H', 'P', 'T');


%% 6) Group difference test: network-level
disp(['## Group difference test: network-level - processing', newline]);
load([outpath, 'wholesub_NET_dc.mat']);
num_network = 7;

H = zeros(num_network, Nstep);
P = zeros(num_network, Nstep);
T = zeros(num_network, Nstep);
for step = 1 : Nstep
    for nidx = 1 : num_network
        [~,p,~,stats] = ttest2(net_dc(group==2,nidx,step), net_dc(group==1,nidx,step));
        P(nidx, step) = p;
        T(nidx, step) = stats.tstat;
    end
    [selected, ~, ~, corrected_P] = fdr_bh(P(:, step), 0.05);
    H(:, step) = selected;
    P(:, step) = corrected_P;
end
save([outpath, 'high_low/groupdiff_NET_ttest.mat'], 'H', 'P', 'T');


%% 7) Group difference test: subcortex-wise
disp(['## Group difference test: subcortex-wise - processing', newline]);
Nsubcor = 7;    % amygdala, hippocampus, globus pallidus, nucleus accumbens, putamen, caudate, and thalamus
subcor_meanDC = zeros(Nsub, Nsubcor, Nstep);
for step = 1 : Nstep
    subcor_dc = mean_subcortical(roi_dc(:,:,step)')';
    % take average for left and right hemisphere
    subcor_meanDC(:,:,step) = mean(reshape(subcor_dc, [Nsub, Nsubcor, 2]), 3);
end

H = zeros(Nsubcor, Nstep);
P = zeros(Nsubcor, Nstep);
T = zeros(Nsubcor, Nstep);
for step = 1 : Nstep
    for ridx = 1 : Nsubcor
        [~,p,~,stats] = ttest2(subcor_meanDC(group==2,ridx,step), subcor_meanDC(group==1,ridx,step));
        P(ridx, step) = p;
        T(ridx, step) = stats.tstat;
    end
    % correct across regions
    [h, ~, ~, p] = fdr_bh(P(:,step), 0.05);
    P(:, step) = p;
    H(:, step) = h;
end
save([outpath, 'high_low/groupdiff_SUB_ttest.mat'], 'H', 'P', 'T');
end
end
