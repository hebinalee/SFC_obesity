%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Perform group analysis using degree of SFC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function group_analysis()
inpath = 'X:/path/myfolder/inputs/';
outpath = 'X:/path/myfolder/outputs/';
Nsub = 301;
Nroi = 246;
Nstep = 5;

%% 1) Divide participants into two group based on obesity phenotypes
load([inpath, 'a_dataset.mat'])
group = obesity(:,1)>=18.5 & obesity(:,1)<25 & ((sex==1 & obesity(:,2)<0.95) | (sex==2 & obesity(:,2)<0.8));    % healthy weight
group = group + 2*obesity(:,1)>=25 & ((sex==1 & obesity(:,2)>1.0) | (sex==2 & obesity(:,2)>0.86));              % overweight
save([outpath, 'a_group.mat'], 'group')


%% 2) Compute group average SFC
load([inpath, 'c_seed_regions.mat']);

grpmean_SFC = cell(2,1);
grpmean_SFC{1} = zeros(length(seed_idx), Nroi, Nstep);
grpmean_SFC{2} = zeros(length(seed_idx), Nroi, Nstep);

for sidx = 1 : Nsub
    disp(['subject = ', num2str(sidx)])
    if group(sidx) == 1
        load([inpath, '2.sfc5/sub', pad(num2str(sidx, '%d'), 3, 'left', '0'), '.mat']);
        grpmean_SFC{1} = grpmean_SFC{1} + sfc(seed_idx,:,:);
    elseif group(sidx) == 2
        load([inpath, '2.sfc5/sub', pad(num2str(sidx, '%d'), 3, 'left', '0'), '.mat']);
        grpmean_SFC{2} = grpmean_SFC{2} + sfc(seed_idx,:,:);
    end
end
grpmean_SFC{1} = grpmean_SFC{1} / sum(group == 1);
grpmean_SFC{2} = grpmean_SFC{2} / sum(group == 2);

grpmean_DC = zeros(2, Nroi, Nstep);
grpmean_DC(1,:,:) = squeeze(sum(grpmean_SFC{1}(:,:,1:Nstep), 1));
grpmean_DC(2,:,:) = squeeze(sum(grpmean_SFC{2}(:,:,1:Nstep), 1));
save([outpath, 'groupmean_SFC.mat'], 'grpmean_SFC', 'grpmean_DC')

    
%% 3) Compute hub regions
hub = cell(2,1);
for gidx = 1 : 2
    hub{gidx} = zeros (224, Nstep);
    for step = 1 : Nstep
        stepDC = squeeze(grpmean_DC(gidx, :, step))';
        stepDC = (stepDC - min(stepDC)) / (max(stepDC) - min(stepDC));
        stepDC = stepDC / mean(stepDC);
        hub{gidx}(:, step) = [stepDC(1:210); mean_subcortical(stepDC)] > 1.5;
    end
end
save([outpath, 'groupdiff_hub.mat'], 'hub');


%% 4) Group difference test: roi-level
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
    % correct across networks
    [selected, ~, ~, corrected_P] = fdr_bh(P(:, step), 0.05);
    H(:, step) = selected;
    P(:, step) = corrected_P;
end
save([outpath, 'groupdiff_ROI_ttest.mat'], 'H', 'P', 'T');


%% 5) Group difference test: network-level
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
save([outpath, 'groupdiff_NET_ttest.mat'], 'H', 'P', 'T');


%% 6) Group difference test: subcortex-wise
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
save([outpath, 'groupdiff_SUB_ttest.mat'], 'H', 'P', 'T');
end
