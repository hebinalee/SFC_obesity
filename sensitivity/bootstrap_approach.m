function bootstrap_approach()
inpath = 'X:/path/myfolder/inputs/';
outpath = 'X:/path/myfolder/outputs/';
Nsub = 301;
Nroi = 246;
Nstep = 7;

%% load DC values & WHR
load([inpath, 'a_dataset.mat'])
load([inpath, 'a_conn_ridge.mat'])
DC = sum(CONN,3);
clear CONN
Nperm = 1000;


%% 1) Subject bootstrapping
subsample = randi(Nsub, round(Nsub*0.9), Nperm);
save([outpath, 'bootstrap/bootstrap_subsamples.mat'], 'subsample')

    
%% 2) Associate degree centrality with WHR
bootstrap_selected = zeros(Nroi, Nperm);
roi_prob_map = zeros(Nroi, 1);
bootstrap_meanR = zeros(Nroi, 1);
bootstrap_meansigR = zeros(Nroi, 1);
P = zeros(Nroi, 1);
for i = 1 : Nperm
    if rem(i,100) == 0
        disp(['i = ', num2str(i)]);
    end
    subX = DC(subsample(:, i), :);      % DC
    subY = obesity(subsample(:, i), 2); % WHR
    for roi = 1 : Nroi
        [r, p] = corrcoef(subX(:, roi), subY);
        P(roi) = p(1,2);
        bootstrap_meanR(roi) = bootstrap_meanR(roi) + r(1,2);
        if bootstrap_selected(roi, i) == 1
            bootstrap_meansigR(roi) = bootstrap_meansigR(roi) + r(1,2);
        end
    end
    [selected, ~, ~, corrected_P] = fdr_bh(P, 0.05);
    bootstrap_selected(:, i) = selected;
    roi_prob_map = roi_prob_map + selected;
end
bootstrap_meanR = bootstrap_meanR / 1000;

for roi = 1 : Nroi
    bootstrap_meansigR(roi) = bootstrap_meansigR(roi) / sum(bootstrap_selected(roi, :));
end
save([outpath, 'bootstrap/bootstrap_seedROI.mat'], 'bootstrap_selected', 'roi_prob_map', 'bootstrap_meanR', 'bootstrap_meansigR')

    
%% 3) Save ROI probability map image
bna = load_nii([store7, 'hblee/ETC/Atlas/BNA/BNA_3mm.nii']);
roi = bna.img;
roi_probmap_img = zeros(size(roi));
meanR_img = zeros(size(roi));
meansigR_img = zeros(size(roi));
for i = 1 : Nroi
    idx = find(roi(:) == i);
    roi_probmap_img(idx) = roi_prob_map(i)/1000;
    meanR_img(idx) = bootstrap_meanR(i);
    meansigR_img(idx) = bootstrap_meansigR(i);
end
bna.img = roi_probmap_img;
save_nii(bna, [outpath, 'bootstrap/bootstrap_seed_probmap.nii']);
bna.img = meanR_img;
save_nii(bna, [outpath, 'bootstrap/bootstrap_meanR.nii']);
bna.img = meansigR_img;
save_nii(bna, [outpath, 'bootstrap/bootstrap_meansigR.nii']);


%% 4) Compute group average SFC matrix
clear DC
Nsubsample = size(subsample, 1);
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
    save([outpath, 'bootstrap/sfc/stepwise_', num2str(i), '.mat'], 'grpmean_SFC')
end


%% 5) average SFC analysis across iteration
bootstrap_mean_sfc = cell(2,1);
bootstrap_mean_sfc{1} = zeros(Nroi, Nstep);
bootstrap_mean_sfc{2} = zeros(Nroi, Nstep);
for i = 1 : Nperm
    load([outpath, 'bootstrap/sfc/stepwise_', num2str(i), '.mat'])
    bootstrap_mean_sfc{1} = bootstrap_mean_sfc{1} + squeeze(sum(seed_avg{1}, 1));
    bootstrap_mean_sfc{2} = bootstrap_mean_sfc{2} + squeeze(sum(seed_avg{2}, 1));
end
bootstrap_mean_sfc{1} = bootstrap_mean_sfc{1} / 1000;
bootstrap_mean_sfc{2} = bootstrap_mean_sfc{2} / 1000;
save([outpath, 'bootstrap/bootstrap_mean_SFC.mat'], 'bootstrap_mean_sfc')
end
