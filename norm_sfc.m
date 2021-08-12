function norm_sfc()
if ispc
    store7 = 'V:/';
else
    store7 = '/store7/';
end

inpath = [store7, 'hblee/2.stepwise/9.regress-301/'];
outpath = [store7, 'hblee/2.stepwise/10.final/'];
Nsub = 301;
Nroi = 246;
Nstep = 200;

%% Subject selection
% group = xlsread([store7, 'hblee/2.stepwise/Enhanced_NKI.xlsx'], 9, 'Q3:Q303');
% save([outpath, 'a_group.mat'], 'group')
load([inpath, 'a_dataset.mat'])
load([outpath, 'a_group.mat'])

%% Regress out from BMI & WHR
% make model
% tbl = table(obesity(group>0, 1), age(group>0), sex(group>0), 'VariableNames', {'bmi','age','sex'});
% mdl = fitlm(tbl, 'bmi ~ age+sex');
% bmi_reg = mdl.Residuals.Raw;
% tbl = table(obesity(group>0, 2), age(group>0), sex(group>0), 'VariableNames', {'whr','age','sex'});
% mdl = fitlm(tbl, 'whr ~ age+sex');
% whr_reg = mdl.Residuals.Raw;
% save([outpath, 'a_obesity_reg.mat'], 'bmi_reg', 'whr_reg')

%% Regress out age & sex & headmotion from CONNECTIVITY
% load([outpath, 'a_dataset.mat'])
% load([inpath, 'a_headmotion.mat'])
% load([inpath, 'a_conn_ridge.mat'])
% CONN_reg = zeros(sum(group>0), Nroi, Nroi);
% for i = 1 : Nroi
%     for j = 1 : Nroi
%         if i ~= j
%             % make model
%             tbl = table(CONN(group>0,i,j), age(group>0), sex(group>0), headmotion(group>0), 'VariableNames', {'R','age','sex', 'headmotion'});
%             mdl = fitlm(tbl, 'R ~ age+sex+headmotion');
%             CONN_reg(:,i,j) = mdl.Residuals.Raw;
%         end
%     end
% end
% save([outpath, 'b_conn_regressed3.mat'], 'CONN_reg', '-v7.3')

%% Compute group average SFC
seed_idx = [7, 8, 9, 20, 27, 28, 54, 58, 59, 60, 89, 93, 99, 109, 152, 162, 165, 198, 212, 222, 231, 233, 234, 237, 238, 243]';
seed_sub{1} = [7, 8, 20, 27, 28, 89, 93, 99, 109, 152, 198]';   % Cognition
seed_sub{2} = [165, 237, 238, 243]';                            % Reward
seed_sub{3} = [9, 54, 58, 59, 60, 162, 222, 231, 233, 234]';    % Sensorimotor

seed_avg = cell(2,1);
seed_avg{1} = zeros(length(seed_idx), Nroi, Nstep);
seed_avg{2} = zeros(length(seed_idx), Nroi, Nstep);
for sidx = 1 : Nsub
    if group(sidx) == 1
        load([inpath, '2.sfc5/sub', pad(num2str(sidx, '%d'), 3, 'left', '0'), '.mat']);
        seed_avg{1} = seed_avg{1} + normalize_sfc(sfc(seed_idx,:,:));
    elseif group(sidx) == 2
        load([inpath, '2.sfc5/sub', pad(num2str(sidx, '%d'), 3, 'left', '0'), '.mat']);
        seed_avg{2} = seed_avg{2} + normalize_sfc(sfc(seed_idx,:,:));
    end
end
seed_avg{1} = seed_avg{1} / sum(group == 1);
seed_avg{2} = seed_avg{2} / sum(group == 2);
save([outpath, 'c_seed_norm_avg.mat'], 'seed_avg')

%% Figure 3) Stepwise connectivity matrix
group_name = {'HW', 'OB'};
for gidx = 1 : 2
    for step = [1, 2, 3, 4, 5, 6, 7] %10, 20, 40, 60, 100, 150, 200]
        imagesc(seed_avg{gidx}(:,:,step), [-1 6]), axis off%, colorbar('FontSize', 24)
%         yticks([1:length(seed_idx)]), yticklabels(seed_idx), xticklabels([])
        saveas(gcf, [outpath, 'figures/norm/sfc/', group_name{gidx}, '-step', num2str(step), '.png']);
    end
end

%% Figure 3) Spider plot of DC
load([inpath, 'cluster_Fan_Net_r280.mat'])
net8 = cluster_Fan_Net.dat(1:246, 3);
net = cluster_Fan_Net.descrip{3,2};
net(9) = [];
steps = [1, 2, 3, 4, 5, 6, 7];
net_dc = zeros(2, length(steps), 8);

load([outpath, 'c_seed_norm_avg.mat'])
for step = steps
    for gidx = 1 : 2
        seedconn = seed_avg{gidx};
        dc = sum(seedconn(:,:,step), 1);
        dc(isinf(dc)|isnan(dc)) = 0;
        for nidx = 1 : 8
            net_dc(gidx, nidx) = mean(dc(net8 == nidx));
        end
    end
    (net_dc(1,:)-net_dc(2,:))/(max(net_dc(:))-min(net_dc(:)))
    % max(a(1,:)-a(2,:))
    spider_plot3(net_dc, 'AxesLabels', net, 'AxesLimits', [repmat(min(net_dc(:)), 1, 8); repmat(max(net_dc(:)), 1, 8)]);
    %if step == 1
    %    legend('Healthy', 'Obesity', 'Location', 'southoutside');
    %end
    saveas(gcf, [outpath, 'figures/norm/network/', num2str(v), '-step', num2str(step),'.png']);
end

%% Figure 3) DC plot
% Surface for 'plot_hemispheres'
parc = zeros(64984,1);
gii = gifti([inpath, 'conte69/BNA_conte69.L.func.gii']);
parc(1:32492,:) = gii.cdata;
gii = gifti([inpath, 'conte69/BNA_conte69.R.func.gii']);
parc(1+32492:end,:) = gii.cdata;
[surf_lh, surf_rh] = load_conte69();

% cortex
load([outpath, 'c_seed_norm_avg.mat'])
groupname = {'HW', 'OB'};
dc_step = zeros(2, Nroi, 4);

steps = [[1, 2, 3, 4]; [5, 6, 7, 8]];
clim = cell(2,1);
% clim = [[-6 -11 -14 -16 -16 -16 -16 -16]; [15 16 22 25 28 28 28 28]];
clim = [repmat(-16, [1 8]); repmat(28, [1 8])];
for mod = 1         % edit here
    for gidx = 1	% edit here
        seedconn = seed_avg{gidx};
        i = 1;
        for step = steps(mod, :)
            dc = sum(seedconn(:,:,step), 1);
            dc(isinf(dc)|isnan(dc)) = 0;
            dc_step(gidx, :,i) = dc;
            i = i + 1;
        end
        obj = plot_hemispheres2(squeeze(dc_step(gidx, 1:210, :)), {surf_lh,surf_rh}, 'parcellation', parc, 'clim', clim(:,4*mod-3:4*mod));
%         saveas(gcf, [outpath, 'figures/norm/dcstep/', groupname{gidx}, num2str(mod), '.png'])
%         close(gcf)
    end
end

% subcortex
% steps = [1, 2, 3, 4, 5, 6, 7];
% for gidx = 1 : 2
%     for step = 1 : length(steps)
%         dc_subcor = mean_subcor(sum(seed_avg{gidx}(:, :, steps(step)))');
%         plot_subcortical(dc_subcor, 'ventricles', 'False', 'cmap', 'parula', 'color_range', clim(:,step)')
%         saveas(gcf, [outpath, 'figures/norm/dcstep/subcor-', groupname{gidx}, num2str(step), '.png'])
% %         close(gcf)
%     end
% end

% colorbar image
% imagesc(seed_avg{1}(:,:,3), [-1 6]), axis off, colorbar('southoutside', 'FontSize', 30)
% imagesc(seed_avg{1}(:,:,3), [-16 28]), axis off
% c = colorbar('southoutside', 'FontSize', 30)
% c.Ticks = [-16,28]

%% Calculate DC vector for all subjects
load([inpath, 'cluster_Fan_Net_r280.mat'])
net8 = cluster_Fan_Net.dat(1:246, 3);
net = cluster_Fan_Net.descrip{3,2};
net(9) = [];

seed_idx = [7, 8, 9, 20, 27, 28, 54, 58, 59, 60, 89, 93, 99, 109, 152, 162, 165, 198, 212, 222, 231, 233, 234, 237, 238, 243]';
seed_sub{1} = [1:26]';
seed_sub{2} = [1, 2, 4, 5, 6, 11, 12, 13, 14, 15, 18]'; % Cognition
seed_sub{3} = [17, 24, 25, 26]';                        % Reward
seed_sub{4} = [3, 7, 8, 9, 10, 16, 20, 21, 22, 23]';    % Sensorimotor

load([outpath, 'a_group.mat'])
steps = [1, 2, 3, 4, 5, 6, 7];
net_dc = zeros(4, Nsub, length(steps), 8);
clear dc
for sidx = 1 : Nsub
    disp(['subject = ', num2str(sidx)])
    load([inpath, '2.sfc5/sub', pad(num2str(sidx, '%d'), 3, 'left', '0'), '.mat'])
    sfc = sfc(seed_idx,:,steps);
    for i = 1 : 4   % seed_type (whole, cog, rew, sen)
        for step = 1 : length(steps)
            dc{i} = sum(normalize_sfc(sfc(seed_sub{i},:,step)), 1);
            dc{i}(isinf(dc{i})|isnan(dc{i})) = 0;
            for nidx = 1 : 8
                net_dc(i, sidx, step, nidx) = mean(dc{i}(net8 == nidx));
            end
        end
    end
end
save([inpath, 'd_wholesub_norm_netdc.mat'], 'net_dc');

HW_dc = net_dc(:, group==1, :, :);
OB_dc = net_dc(:, group==2, :, :);
save([outpath, 'd_group_norm_netdc.mat'], 'HW_dc', 'OB_dc');

%% Group difference - network level
H = zeros(length(steps), 8, 1);
P = zeros(length(steps), 8, 1);
T = zeros(length(steps), 8, 1);
for i = 1   % seed_type (whole, cog, rew, sen)
    for step = 1 : length(steps)
        for nidx = 1 : 8
            [~,p,~,stats] = ttest2(squeeze(OB_dc(i,:,step,nidx)),squeeze(HW_dc(i,:,step,nidx)));
%             H(step,nidx,i) = h;
            P(step,nidx,i) = p;
            T(step,nidx,i) = stats.tstat;
        end
    end
end
[H, ~, ~, P] = fdr_bh(P, 0.05);
save([outpath, 'e_net_ttest.mat'], 'H', 'P', 'T');

%% Group difference - ROI level
seed_idx = [7, 8, 9, 20, 27, 28, 54, 58, 59, 60, 89, 93, 99, 109, 152, 162, 165, 198, 212, 222, 231, 233, 234, 237, 238, 243]';
seed_sub{1} = [1:26]';
seed_sub{2} = [1, 2, 4, 5, 6, 11, 12, 13, 14, 15, 18]'; % Cognition
seed_sub{3} = [17, 24, 25, 26]';                        % Reward
seed_sub{4} = [3, 7, 8, 9, 10, 16, 20, 21, 22, 23]';    % Sensorimotor

load([outpath, 'a_group.mat'])
steps = [1, 2, 3, 4, 5, 6, 7];
net_dc = zeros(4, Nsub, Nroi, length(steps));
clear dc
for sidx = 1 : Nsub
    disp(['subject = ', num2str(sidx)])
    load([inpath, '2.sfc5/sub', pad(num2str(sidx, '%d'), 3, 'left', '0'), '.mat'])
    sfc = sfc(seed_idx,:,steps);
    for i = 1 : 4   % seed_type (whole, cog, rew, sen)
        for step = 1 : length(steps)
            dc{i} = sum(normalize_sfc(sfc(seed_sub{i},:,step)), 1);
            dc{i}(isinf(dc{i})|isnan(dc{i})) = 0;
            net_dc(i, sidx, :, step) = dc{i};
        end
    end
end
save([inpath, 'd_wholesub_norm_roidc.mat'], 'net_dc');

HW_dc = net_dc(:, group==1, :, :);
OB_dc = net_dc(:, group==2, :, :);
save([outpath, 'd_group_norm_roidc.mat'], 'HW_dc', 'OB_dc');
%%
Nroi = 246;
H = zeros(4, Nroi, length(steps));
P = zeros(4, Nroi, length(steps));
T = zeros(4, Nroi, length(steps));

for i = 1 : 4
    for step = 1 : length(steps)
%         ob_dc = [squeeze(OB_dc(i,:,1:210,1)), mean_subcor(squeeze(OB_dc(i,:,:,1))')'];
%         hw_dc = [squeeze(HW_dc(i,:,1:210,1)), mean_subcor(squeeze(HW_dc(i,:,:,1))')'];
        for roi = 1 : Nroi
            [~,p,~,stats] = ttest2(squeeze(OB_dc(i,:,roi,step)), squeeze(HW_dc(i,:,roi,step)));
%             [~,p,~,stats] = ttest2(squeeze(ob_dc(:,roi)), squeeze(hw_dc(:,roi)));
            P(i, roi, step) = p;
            T(i, roi, step) = stats.tstat;
        end
        [selected, ~, ~, corrected_P] = fdr_bh(P(i, :, step), 0.05);
        H(i, :, step) = selected;
        P(i, :, step) = corrected_P;
    end
end
save([outpath, 'e_roi_ttest_', num2str(Nroi), '.mat'], 'H', 'P', 'T');

%%
% Plot p-values for significant regions
% significantP =  H .* P .* sign(T);
significantP =  H .* T;
load([inpath, 'cool_warm.mat'])

seed_type = 1;  % edit here
for mod = 2     % edit here
    if mod == 1
        pp = squeeze(significantP(seed_type, 1:210, 1:4, 1));
%         clim = [repmat(-0.05,1,4); repmat(0.05,1,4)];
        clim = [repmat(-5,1,4); repmat(5,1,4)];
    else
        pp = squeeze(significantP(seed_type, 1:210, 5:7, 1));
%         clim = [repmat(-0.05,1,3); repmat(0.05,1,3)];
        clim = [repmat(-5,1,3); repmat(5,1,3)];
    end
    obj = plot_hemispheres2(pp, {surf_lh,surf_rh}, 'parcellation', parc, 'clim', clim);
%     colormap(obj.figure, enigma_colormap('RdBu_r'));
    colormap(obj.figure, [cool_warm(129,:); cool_warm]);
%     saveas(gcf, [outpath, 'figures/norm/dcstep/', groupname{gidx}, num2str(mod), '.png'])
%     close(gcf)
end
%%
steps = [1, 2, 3, 4, 5, 6, 7];
for i = 3 : 4
%     p_subcor = mean_subcor(squeeze(significantP(i, :, :)));
    p_subcor = squeeze(significantP(i, 211:end, :));
    for step = 1 : length(steps)
        plot_subcortical(p_subcor(:,step), 'ventricles', 'False', 'cmap', 'RdBu_r', 'color_range', [-5 5])
        enigma_colormap(cool_warm)
        saveas(gcf, [outpath, 'figures/norm/pvalue/buf/subcor-T', num2str(i), '-', num2str(step), '.png'])
%         close(gcf)
    end
end

% colorbar image
% imagesc(seed_avg{1}(:,:,3), [-0.05 0.05]), axis off, colormap(cool_warm)
% c = colorbar('southoutside', 'FontSize', 30)
% c.Ticks = [-0.05, 0.05]
 
%% S.Figure 1) Cognition-related seed
clim{1} = [[0 1 14 0.2e4 0.03e11 0.08e23 0.01e49 0.0001e76]; [3 30 480 14e4 2.3e11 11e23 15e49 3.4e76]];
clim{2} = [[0 1 14 0.2e4 0.03e11 0.08e23 0.1e49 0.1e74]; [3 30 480 14e4 2.3e11 8.5e23 3e49 17e74]];

load([outpath, 'cog/c_seed_avg.mat'])
groupname = {'HW', 'OB'};
dc_step = zeros(2, Nroi, 4);
% cortex
% mod = 2;	% edit here
% steps = [[1, 2, 3, 5]; [10, 20, 40, 60]];
% for gidx = 1	% edit here
%     seedconn = seed_avg{gidx};    
%     i = 1;
%     for step = steps(mod, :)
%         dc = sum(seedconn(:,:,step), 1);
%         dc(isinf(dc)|isnan(dc)) = 0;
%         dc_step(gidx, :,i) = dc;
%         i = i + 1;
%     end
%     obj = plot_hemispheres2(squeeze(dc_step(gidx, 1:210, :)), {surf_lh,surf_rh}, 'parcellation', parc, 'clim', clim{gidx}(:,4*mod-3:4*mod));
% end

% subcortex
% steps = [1, 2, 3, 5, 10, 20, 40, 60];
% for gidx = 1 : 2
%     for step = 1 : 8
%         dc_subcor = mean_subcor(sum(seed_avg{gidx}(:, :, steps(step)))');
%         plot_subcortical(dc_subcor, 'ventricles', 'False', 'cmap', 'parula', 'color_range', clim{gidx}(:,step)')
%         saveas(gcf, [outpath, 'figure/paper/supplement/sfig1-cog-subcor-', groupname{gidx}, num2str(step), '.png'])
%         close(gcf)
%     end
% end

% Stepwise connectivity matrix
for gidx = 1 : 2
    for step = [1, 2, 3, 5, 10, 20, 40, 60]
        imagesc(seed_avg{gidx}(:,:,step)), axis off, colorbar('FontSize', 20)
%         yticks([1:length(seed_idx)]), yticklabels(seed_idx), xticklabels([])
        saveas(gcf, [outpath, 'figure/sfc/supplement/cog-', group_name{gidx}, '-step', num2str(step), '.png']);
    end
end

%% S.Figure 2) Reward-related seed
clim{1} = [[0 0.02 2 490 0.03e10 0.1e22 0.01e48 0.0002e75]; [1.1 3.5 50 13400 2.3e10 12.5e22 19e48 4.2e75]];
clim{2} = [[0 0.02 2 490 0.03e10 0.1e22 0.09e48 0.2e73]; [1.1 3.5 50 13400 2.3e10 8.5e22 2.5e48 11.3e73]];

load([outpath, 'rew/c_seed_avg.mat'])
groupname = {'HW', 'OB'};
dc_step = zeros(2, Nroi, 4);
% cortex
% mod = 2;	% edit here
% steps = [[1, 2, 3, 5]; [10, 20, 40, 60]];
% for gidx = 1	% edit here
%     seedconn = seed_avg{gidx};    
%     i = 1;
%     for step = steps(mod, :)
%         dc = sum(seedconn(:,:,step), 1);
%         dc(isinf(dc)|isnan(dc)) = 0;
%         dc_step(gidx, :,i) = dc;
%         i = i + 1;
%     end
%     obj = plot_hemispheres2(squeeze(dc_step(gidx, 1:210, :)), {surf_lh,surf_rh}, 'parcellation', parc, 'clim', clim{gidx}(:,4*mod-3:4*mod));
% end

% subcortex
% steps = [1, 2, 3, 5, 10, 20, 40, 60];
% for gidx = 1 : 2
%     for step = 1 : 8
%         dc_subcor = mean_subcor(sum(seed_avg{gidx}(:, :, steps(step)))');
%         plot_subcortical(dc_subcor, 'ventricles', 'False', 'cmap', 'parula', 'color_range', clim{gidx}(:,step)')
%         saveas(gcf, [outpath, 'figure/paper/supplement/sfig2-rew-subcor-', groupname{gidx}, num2str(step), '.png'])
%         close(gcf)
%     end
% end

% Stepwise connectivity matrix
for gidx = 1 : 2
    for step = [1, 2, 3, 5, 10, 20, 40, 60]
        imagesc(seed_avg{gidx}(:,:,step)), axis off, colorbar('FontSize', 20)
%         yticks([1:length(seed_idx)]), yticklabels(seed_idx), xticklabels([])
        saveas(gcf, [outpath, 'figure/sfc/supplement/rew-', group_name{gidx}, '-step', num2str(step), '.png']);
    end
end

%% S.Figure 3) Sensorimotor-related seed
clim{1} = [[0 0.1 3 0.09e4 0.1e10 0.03e23 0.005e49 0.0004e75]; [3.6 26 300 4.6e4 8.3e10 5.3e23 9e49 20e75]];
clim{2} = [[0 0.1 3 0.09e4 0.1e10 0.03e23 0.1e48 0.04e74]; [3.6 26 300 4.6e4 8.3e10 3e23 11e48 6.5e74]];

load([outpath, 'sen/c_seed_avg.mat'])
groupname = {'HW', 'OB'};
dc_step = zeros(2, Nroi, 4);
% cortex
% mod = 2;	% edit here
% steps = [[1, 2, 3, 5]; [10, 20, 40, 60]];
% for gidx = 2	% edit here
%     seedconn = seed_avg{gidx};    
%     i = 1;
%     for step = steps(mod, :)
%         dc = sum(seedconn(:,:,step), 1);
%         dc(isinf(dc)|isnan(dc)) = 0;
%         dc_step(gidx, :,i) = dc;
%         i = i + 1;
%     end
%     obj = plot_hemispheres2(squeeze(dc_step(gidx, 1:210, :)), {surf_lh,surf_rh}, 'parcellation', parc, 'clim', clim{gidx}(:,4*mod-3:4*mod));
% end

% subcortex
% steps = [1, 2, 3, 5, 10, 20, 40, 60];
% for gidx = 1 : 2
%     for step = 1 : 8
%         dc_subcor = mean_subcor(sum(seed_avg{gidx}(:, :, steps(step)))');
%         plot_subcortical(dc_subcor, 'ventricles', 'False', 'cmap', 'parula', 'color_range', clim{gidx}(:,step)')
%         saveas(gcf, [outpath, 'figure/paper/supplement/sfig3-sen-subcor-', groupname{gidx}, num2str(step), '.png'])
%         close(gcf)
%     end
% end

% Stepwise connectivity matrix
for gidx = 1 : 2
    for step = [1, 2, 3, 5, 10, 20, 40, 60]
        imagesc(seed_avg{gidx}(:,:,step)), axis off, colorbar('FontSize', 20)
%         yticks([1:length(seed_idx)]), yticklabels(seed_idx), xticklabels([])
        saveas(gcf, [outpath, 'figure/sfc/supplement/sen-', group_name{gidx}, '-step', num2str(step), '.png']);
    end
end
end