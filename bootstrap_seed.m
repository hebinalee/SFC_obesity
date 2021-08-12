function bootstrap_seed()
if ispc
    store7 = 'V:/';
else
    store7 = '/store7/';
end
inpath = [store7, 'hblee/2.stepwise/9.regress-301/'];
outpath = [store7, 'hblee/2.stepwise/10.final/'];

%% load DC values & WHR
load([inpath, 'a_dataset.mat'])
load([inpath, 'a_conn_ridge.mat'])
DC = sum(CONN,3);
clear CONN
Nsub = size(DC, 1);
Nroi = size(DC, 2);
Nperm = 1000;

%% Subject bootstrapping
subsample = randi(Nsub, round(Nsub*0.9), Nperm);
save([outpath, 'revision/bootstrap_subsamples.mat'], 'subsample')

%% Calculate ROI probability map
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
save([outpath, 'revision/bootstrap_seedROI.mat'], 'bootstrap_selected', 'roi_prob_map', 'bootstrap_meanR', 'bootstrap_meansigR')

%% Save ROI probability map image
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
save_nii(bna, [outpath, 'revision/bootstrap_seed_probmap.nii']);
bna.img = meanR_img;
save_nii(bna, [outpath, 'revision/bootstrap_meanR.nii']);
bna.img = meansigR_img;
save_nii(bna, [outpath, 'revision/bootstrap_meansigR.nii']);

%% Report subcortical results
% Surface for 'plot_hemispheres'
parc = zeros(64984,1);
gii = gifti([inpath, 'conte69/BNA_conte69.L.func.gii']);
parc(1:32492,:) = gii.cdata;
gii = gifti([inpath, 'conte69/BNA_conte69.R.func.gii']);
parc(1+32492:end,:) = gii.cdata;
[surf_lh, surf_rh] = load_conte69();

load([inpath, 'cool_warm.mat'])
cool_warm(1,:) = cool_warm(129,:);

% cortex
obj = plot_hemispheres2(bootstrap_meanR(1:210,:), {surf_lh,surf_rh}, 'parcellation', parc, 'clim', [-0.3; 0.3]);
colormap(obj.figure, cool_warm)
saveas(gcf, [outpath, 'revision/figures/bootstrap_R.png'])
close(gcf)

% subcortical
subcor_meanR = mean_subcor(bootstrap_meanR);
plot_subcortical(subcor_meanR, 'ventricles', 'False', 'cmap', 'RdBu_r', 'color_range', [-0.3 0.3]);
enigma_colormap(cool_warm)
saveas(gcf, [outpath, 'revision/figures/bootstrap_R_subcor.png'])
close(gcf)

%% SFC analysis
% clear DC
Nstep = 7;
Nsubsample = size(subsample, 1);
for i = 1 : Nperm
    if rem(i,100) == 0
        disp(['i = ', num2str(i)]);
    end
    seed_idx = find(bootstrap_selected(:,i)==1);    
    seed_avg = cell(2,1);
    seed_avg{1} = zeros(length(seed_idx), Nroi, Nstep);
    seed_avg{2} = zeros(length(seed_idx), Nroi, Nstep);
    for j = 1 : Nsubsample
        if group(subsample(j,i)) == 1
            load([inpath, '2.sfc5/sub', pad(num2str(subsample(j,i), '%d'), 3, 'left', '0'), '.mat']);
            seed_avg{1} = seed_avg{1} + normalize_sfc(sfc(seed_idx,:,1:Nstep));
        elseif group(subsample(j,i)) == 2
            load([inpath, '2.sfc5/sub', pad(num2str(subsample(j,i), '%d'), 3, 'left', '0'), '.mat']);
            seed_avg{2} = seed_avg{2} + normalize_sfc(sfc(seed_idx,:,1:Nstep));
        end
    end
    seed_avg{1} = seed_avg{1} / sum(group(subsample(:,i)) == 1);
    seed_avg{2} = seed_avg{2} / sum(group(subsample(:,i)) == 2);
    save([outpath, 'revision/bootstrap_sfc/stepwise_', num2str(i), '.mat'], 'seed_avg')
end

%% average SFC analysis
bootstrap_mean_sfc = cell(2,1);
bootstrap_mean_sfc{1} = zeros(Nroi, Nstep);
bootstrap_mean_sfc{2} = zeros(Nroi, Nstep);
for i = 1 : Nperm
    load([outpath, 'revision/bootstrap_sfc/stepwise_', num2str(i), '.mat'])
    bootstrap_mean_sfc{1} = bootstrap_mean_sfc{1} + squeeze(sum(seed_avg{1}, 1));
    bootstrap_mean_sfc{2} = bootstrap_mean_sfc{2} + squeeze(sum(seed_avg{2}, 1));
end
bootstrap_mean_sfc{1} = bootstrap_mean_sfc{1} / 1000;
bootstrap_mean_sfc{2} = bootstrap_mean_sfc{2} / 1000;
save([outpath, 'revision/bootstrap_mean_SFC.mat'], 'bootstrap_mean_sfc')

%% Plot bootstrap SFC results on hemisphere
parc = zeros(64984,1);
gii = gifti([inpath, 'conte69/BNA_conte69.L.func.gii']);
parc(1:32492,:) = gii.cdata;
gii = gifti([inpath, 'conte69/BNA_conte69.R.func.gii']);
parc(1+32492:end,:) = gii.cdata;
[surf_lh, surf_rh] = load_conte69();

% cortex
% load([outpath, 'revision/bootstrap_mean_SFC.mat'])
group_name = {'HW', 'OB'};
dc_step = zeros(2, Nroi, 4);

steps = [[1, 2, 3, 4]; [5, 6, 7, 7]];
% clim = [[-7 -18 -23 -26 -28 -28 -28 -28]; [10 20 30 40 47 47 47 47]];
clim = [[-10 -20 -30 -30 -30 -30 -30 -30]; [10 20 30 40 50 50 50 50]];
load([inpath, 'cool_warm.mat'])
cool_warm(1,:)=cool_warm(129,:);

for mod = 2
         % edit here
    for gidx = 1	% edit here
        seedconn = bootstrap_mean_sfc{gidx};
        i = 1;
        for step = steps(mod, :)
            dc = seedconn(:,step);
            dc(isinf(dc)|isnan(dc)) = 0;
            dc_step(gidx, :,i) = dc;
            i = i + 1;
        end
        obj = plot_hemispheres2(squeeze(dc_step(gidx, 1:210, :)), {surf_lh,surf_rh}, 'parcellation', parc, 'clim', clim(:,4*mod-3:4*mod));
        colormap(obj.figure, cool_warm)
%         obj.figure.Colormap = cool_warm;
        saveas(gcf, [outpath, 'revision/figures/bootstrap_sfc/', group_name{gidx}, num2str(mod), '.png'])
%         close(gcf)
    end
end
% subcortex
Nstep = 7;
for gidx = 1 : 2
    for step = 1 : Nstep
        dc_subcor = mean_subcor(bootstrap_mean_sfc{gidx}(:, step));
        plot_subcortical(dc_subcor, 'ventricles', 'False', 'cmap', 'RdBu_r', 'color_range', clim(:,step)')
        enigma_colormap(cool_warm)
        saveas(gcf, [outpath, 'revision/figures/bootstrap_sfc/subcor_', group_name{gidx}, num2str(step), '.png'])
%         close(gcf)
    end
end

%% Hub regions
hub = cell(2,1);
hub_HW = bootstrap_mean_sfc{1};
hub_HW = hub_HW ./ mean(hub_HW,1);   % hub_HW = bsxfun(@rdivide, hub_HW, mean(hub_HW));
hub{1} = [hub_HW(1:210, :); mean_subcor(hub_HW)] > 1.5;
hub_OB = bootstrap_mean_sfc{2};
hub_OB = hub_OB ./ mean(hub_OB,1);
hub{2} = [hub_OB(1:210, :); mean_subcor(hub_OB)] > 1.5;
clear hub_HW hub_OB

% cortex
% clim = [[-2 2]; [-2 2]; [-2 2]; [-2 2]];
for gidx = 2
    for step = 2
        if step == 1
            hub_4steps = hub{gidx}(1:210, 1:4);
        else
            hub_4steps = hub{gidx}(1:210, 5:7);
        end
        % remove colormap from plot_hemispheres2
        obj = plot_hemispheres2(hub_4steps, {surf_lh,surf_rh}, 'parcellation', parc);%, 'clim', [-3 3]);
        colormap(obj.figure, cool_warm(129:230, :))
%         obj.figure.Colormap = cool_warm;
        saveas(gcf, [outpath, 'revision/figures/bootstrap_sfc/hub_', group_name{gidx}, num2str(step), '.png'])
        close(gcf)
    end
end
% subcortex
for gidx = 1 : 2
    for step = 1 : 7
        hub_1step = hub{gidx}(:, step);
        plot_subcortical(hub_1step(211:end), 'ventricles', 'False', 'color_range', [-1.7 1.7])
        enigma_colormap(cool_warm)
        saveas(gcf, [outpath, 'revision/figures/bootstrap_sfc/hub_subcor_', group_name{gidx}, num2str(step), '.png'])
    end
end

%% T-test
H = zeros(Nroi, Nstep);
P = zeros(Nroi, Nstep);
T = zeros(Nroi, Nstep);

for step = 1 : Nstep
%     ob_dc = [bootstrap_mean_sfc{2}(1:210,:), mean_subcor(bootstrap_mean_sfc{1})'];
%     hw_dc = [bootstrap_mean_sfc{1}(1:210,:), mean_subcor(bootstrap_mean_sfc{2})'];
    for roi = 1 : Nroi
        [~,p,~,stats] = ttest2(bootstrap_mean_sfc{1}(roi,step), bootstrap_mean_sfc{2}(roi,step));
%         [~,p,~,stats] = ttest2(squeeze(ob_dc(:,roi)), squeeze(hw_dc(:,roi)));
        P(roi, step) = p;
        T(roi, step) = stats.tstat;
    end
    [selected, ~, ~, corrected_P] = fdr_bh(P(i, :, step), 0.05);
    H(i, :, step) = selected;
    P(i, :, step) = corrected_P;
    save([outpath, 'revision/e_roi_ttest_', num2str(Nroi), '.mat'], 'H', 'P', 'T');
end
end
