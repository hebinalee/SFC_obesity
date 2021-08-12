function high_low_risk()
if ispc
    store7 = 'V:/';
else
    store7 = '/store7/';
end

inpath = [store7, 'hblee/2.stepwise/9.regress-301/'];
outpath = [store7, 'hblee/2.stepwise/10.final/'];
Nsub = 301;
Nroi = 246;
Nstep = 7;

%% Subject selection
load([inpath, 'a_dataset.mat'])
group = (sex==1 & obesity(:,2)<0.95) | (sex==2 & obesity(:,2)<0.8);             % low-risk
group = group + 2*((sex==1 & obesity(:,2)>1.0) | (sex==2 & obesity(:,2)>0.86)); % high-risk
% low : high = 184 : 51

%% Compute group average SFC
seed_idx = [7, 8, 9, 20, 27, 28, 54, 58, 59, 60, 89, 93, 99, 109, 152, 162, 165, 198, 212, 222, 231, 233, 234, 237, 238, 243]';

grpmean_SFC = cell(2,1);
grpmean_SFC{1} = zeros(length(seed_idx), Nroi, Nstep);
grpmean_SFC{2} = zeros(length(seed_idx), Nroi, Nstep);
for sidx = 1 : Nsub
    disp(['subject = ', num2str(sidx)])
    if group(sidx) == 1
        load([inpath, '2.sfc5/sub', pad(num2str(sidx, '%d'), 3, 'left', '0'), '.mat']);
        grpmean_SFC{1} = grpmean_SFC{1} + normalize_sfc(sfc(seed_idx,:,:));
    elseif group(sidx) == 2
        load([inpath, '2.sfc5/sub', pad(num2str(sidx, '%d'), 3, 'left', '0'), '.mat']);
        grpmean_SFC{2} = grpmean_SFC{2} + normalize_sfc(sfc(seed_idx,:,:));
    end
end
grpmean_SFC{1} = grpmean_SFC{1} / sum(group == 1);
grpmean_SFC{2} = grpmean_SFC{2} / sum(group == 2);

grpmean_DC = zeros(2, Nroi, Nstep);
grpmean_DC(1,:,:) = squeeze(sum(grpmean_SFC{1}(:,:,1:Nstep), 1));
grpmean_DC(2,:,:) = squeeze(sum(grpmean_SFC{2}(:,:,1:Nstep), 1));
save([outpath, 'revision/highlow_groupmeanSFC.mat'], 'grpmean_SFC', 'grpmean_DC')

%% 1) DC plot
% Surface for 'plot_hemispheres'
parc = zeros(64984,1);
gii = gifti([inpath, 'conte69/BNA_conte69.L.func.gii']);
parc(1:32492,:) = gii.cdata;
gii = gifti([inpath, 'conte69/BNA_conte69.R.func.gii']);
parc(1+32492:end,:) = gii.cdata;
[surf_lh, surf_rh] = load_conte69();

load([inpath, 'cool_warm.mat'])
cool_warm(1,:) = cool_warm(129,:);

% load([outpath, 'revision/highlow_groupmeanSFC.mat'])
% clear grpmean_SFC
group_name = {'HW', 'OB'};

% cortex
for i = 1 : Nstep
    dc_range(i,:)=[min(grpmean_DC(1,1:210,i)), max(grpmean_DC(1,1:210,i)), min(grpmean_DC(2,1:210,i)), max(grpmean_DC(2,1:210,i))];
end

steps = [[1, 2, 3, 4]; [5, 6, 7, 7]];
% clim = cell(2,1);
% clim{1} = [[0 2 26 0.5e4 0.08e11 0.03e24 0.02e50 0.06e76]; [4.3 37 660 21e4 3.7e11 2.7e24 5.5e50 18e76]];
% clim{2} = [[0 2 26 0.5e4 0.16e11 0.1e24 0.01e51 0.0002e79]; [4.3 37 660 21e4 4.7e11 6.5e24 13e51 3.4e79]];
for mod = 1 : 2         % edit here
    for gidx = 1 : 2	% edit here
        dc = squeeze(grpmean_DC(gidx,:,steps(mod,:)));
        dc(isinf(dc)|isnan(dc)) = 0;
        obj = plot_hemispheres2(squeeze(dc(1:210,:)), {surf_lh,surf_rh}, 'parcellation', parc);%, 'clim', clim{gidx}(:,4*mod-3:4*mod));
        colormap(obj.figure, cool_warm)
        saveas(gcf, [outpath, 'revision/figures/highlow/dc_', group_name{gidx}, num2str(mod), '.png'])
        close(gcf)
    end
end

%% subcortex
for gidx = 1 : 2
    for step = 1 : Nstep
        dc_subcor = mean_subcor(squeeze(grpmean_DC(gidx, :, step))');
        plot_subcortical(dc_subcor, 'ventricles', 'False', 'cmap', 'RdBu_r', 'color_range', dc_range(step,gidx*2-1:gidx*2));%clim{gidx}(:,step)')
        enigma_colormap(cool_warm)
        saveas(gcf, [outpath, 'revision/figures/highlow/dc_subcor_', group_name{gidx}, num2str(step), '.png'])
%         close(gcf)
    end
end

%% 2) Hub regions
load([outpath, 'revision/highlow_groupmeanSFC.mat'])
% clear grpmean_SFC
hub = cell(2,1);
for gidx = 1 : 2
    hub{gidx} = zeros (224, Nstep);
    for step = 1 : Nstep
        stepDC = squeeze(grpmean_DC(gidx, :, step))';
        stepDC = (stepDC - min(stepDC)) / (max(stepDC) - min(stepDC));
        stepDC = stepDC / mean(stepDC);
        hub{gidx}(:, step) = [stepDC(1:210); mean_subcor(stepDC)] > 1.5;
    end
end

% cortex
clim = [[-2 -2 -2 -2]; [2 2 2 2]];
for gidx = 1 : 2
    for mode = 1 : 2
        if mode == 1
            hub_4steps = hub{gidx}(1:210, 1:4);
        else
            hub_4steps = hub{gidx}(1:210, 5:7);
        end
        % remove colormap from plot_hemispheres2
        obj = plot_hemispheres2(hub_4steps, {surf_lh,surf_rh}, 'parcellation', parc);%, 'clim', [-3 3]);
        colormap(obj.figure, cool_warm(129:230, :))
        saveas(gcf, [outpath, 'revision/figures/highlow/hub_', group_name{gidx}, num2str(mode), '.png'])
        close(gcf)
    end
end

% subcortex  --  only found in step 1
for gidx = 1 : 2
    for step = 1
        hub_1step = hub{gidx}(:, step);
        plot_subcortical(hub_1step(211:end), 'ventricles', 'False', 'cmap', 'RdBu_r', 'color_range', [-1.3 1.3])
        enigma_colormap(cool_warm)
        saveas(gcf, [outpath, 'revision/figures/highlow/hub_subcor_', group_name{gidx}, num2str(step), '.png'])
    end
end

%% 3) Group difference: roi-level
load([outpath, 'd_wholesub_norm_roidc.mat']);
net_dc = squeeze(net_dc(1,:,:,:));
Nroi = 224;
H = zeros(Nroi, Nstep);
P = zeros(Nroi, Nstep);
T = zeros(Nroi, Nstep);

for step = 1 : Nstep
    ob_dc = [net_dc(group==2,1:210,step), mean_subcor(net_dc(group==2,:,step)')'];
    hw_dc = [net_dc(group==1,1:210,step), mean_subcor(net_dc(group==1,:,step)')'];
    for roi = 1 : Nroi
%         [~,p,~,stats] = ttest2(net_dc(group==2,roi,step), net_dc(group==1,roi,step));
        [~,p,~,stats] = ttest2(ob_dc(:,roi), hw_dc(:,roi));
        P(roi, step) = p;
        T(roi, step) = stats.tstat;
    end
    [selected, ~, ~, corrected_P] = fdr_bh(P(:, step), 0.05);
    H(:, step) = selected;
    P(:, step) = corrected_P;
end
save([outpath, 'revision/highlow_roi_ttest.mat'], 'H', 'P', 'T');

% Plot results
sigT = H .* T;

steps = [[1, 2, 3, 4]; [5, 6, 7, 7]];
clim = [[-5 -5 -5 -5]; [5 5 5 5]];
for mod = 1 : 2
    obj = plot_hemispheres2(sigT(1:210,steps(mod,:)), {surf_lh,surf_rh}, 'parcellation', parc, 'clim', clim);
    colormap(obj.figure, cool_warm)
    saveas(gcf, [outpath, 'revision/figures/highlow/ROIttest', num2str(mod), '.png'])
    close(gcf)
end

% subcortex -- no region survived
for gidx = 1 : 2
    for step = 1
        plot_subcortical(sigT(211:end, step), 'ventricles', 'False', 'cmap', 'RdBu_r', 'color_range', [-5 5])
        enigma_colormap(cool_warm)
        saveas(gcf, [outpath, 'revision/figures/highlow/ROIttest_subcor_', group_name{gidx}, num2str(step), '.png'])
%         close(gcf)
    end
end

%% 4) Group difference: network-level
load([outpath, 'd_wholesub_norm_netdc.mat']);
net_dc = squeeze(net_dc(1,:,:,:));
num_network = 8;

H = zeros(Nstep, num_network);
P = zeros(Nstep, num_network);
T = zeros(Nstep, num_network);
for step = 1 : Nstep
    for nidx = 1 : num_network
        [~,p,~,stats] = ttest2(net_dc(group==2,step,nidx), net_dc(group==1,step,nidx));
%         H(step,nidx,i) = h;
        P(step, nidx) = p;
        T(step, nidx) = stats.tstat;
    end
    [selected, ~, ~, corrected_P] = fdr_bh(P(step, :), 0.05);
    H(step, :) = selected;
    P(step, :) = corrected_P;
end

%Plot bar graph
load([inpath, 'cool_warm.mat'])
red = cool_warm(129+128/2, :);
blue = cool_warm(128/2, :);
for nidx = 1 : num_network
    clf, figure(1)
    hold on
    for step = 1 : Nstep
        h = bar(step, T(step, nidx));
%         if H(i, ridx) == 1
        if T(step, nidx) > 0
            set(h, 'FaceColor', red);
        else
            set(h, 'FaceColor', blue);
        end
    end
    set(gca,'FontSize',15), ylim([-4, 4]), yticks([-4:1:4]), xticks([1:1:7])
%     xlabel('Steps'), ylabel('t-statitic')
    hold off
    saveas(gcf, [outpath, 'revision/figures/highlow/NETttest_', num2str(nidx),'.png']);
end
end