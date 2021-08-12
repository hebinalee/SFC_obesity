function control_confound()
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
% tbl = table(obesity(:, 1), age, sex, 'VariableNames', {'bmi','age','sex'});
% mdl = fitlm(tbl, 'bmi ~ age+sex');
% bmi_reg = mdl.Residuals.Raw;
% tbl = table(obesity(:, 2), age, sex, 'VariableNames', {'whr','age','sex'});
% mdl = fitlm(tbl, 'whr ~ age+sex');
% whr_reg = mdl.Residuals.Raw;
% save([outpath, 'revision/control_obesity_reg.mat'], 'bmi_reg', 'whr_reg')

%% 1) Regress out age & sex & headmotion from CONNECTIVITY
% load([outpath, 'a_dataset.mat'])
load([inpath, 'a_headmotion.mat'])
load([inpath, 'a_conn_ridge.mat'])
CONN_reg = zeros(Nsub, Nroi, Nroi);
for i = 1 : Nroi
    disp(['i = ', num2str(i)])
    for j = 1 : Nroi
        if i ~= j
            % make model
            tbl = table(CONN(:,i,j), age, sex, 'VariableNames', {'R','age','sex'});
            mdl = fitlm(tbl, 'R ~ age+sex');
            CONN_reg(:,i,j) = mdl.Residuals.Raw;
        end
    end
end
save([outpath, 'revision/control_CONN.mat'], 'CONN_reg', '-v7.3')

%% 2) Find WHR-correlated ROIs (DEFINE SEED) -- No seed identified
% load([inpath, 'a_dataset.mat'])
% % load([outpath, 'revision/control_conn.mat'])
% DC = sum(CONN_reg, 3);
% for roi = 1 : Nroi
% %     [r, p] = corrcoef(obesity(:,1), DC(:,roi));
% %     [r, p] = corrcoef(whr_reg, DC(:,roi));
%     R(roi, 1) = r(1,2);
%     P(roi, 1) = p(1,2);
% end
% [selected, ~, ~, corrected_P] = fdr_bh(P, 0.05);
% seed_idx = find(selected == 1);
% significant_R = R(seed_idx);
% significant_P = corrected_P(seed_idx);
% R_seed = selected .* R;
% disp(['p1 = ', num2str(sum(P<0.05)), ' & p2 = ', num2str(length(seed_idx))])
% save([outpath, 'revision/control_Rvalue.mat'], 'R', 'R_seed')

%% Visualize seed on surface
bna = load_nii('Z:/hblee/ETC/Atlas/BNA/BNA_3mm.nii');
roi = bna.img;

WHR_related = zeros(size(roi));
for i = 1 : length(seed_idx)
    idx = find(roi(:) == seed_idx(i));
    WHR_related(idx) = significant_R(i);
end

%% 3) Binarize connectivity
disp(['Binarization of connectivity - processing', newline]);
opt_sparsity = zeros(Nsub, 1);
th = 5;
binCONN = zeros(Nsub, Nroi, Nroi);

for sidx = 1 : Nsub
    disp(['subject = ', num2str(sidx)])
    % load([basepath, 'revision/1.conn0.5/conn-sub', pad(num2str(sidx, '%d'), 3, 'left', '0'), '.mat'])
    conn = squeeze(CONN_reg(sidx, :, :));
    binconn = zeros(Nroi, Nroi);
    W = conn;
    if isequal(W,W.')	% if symmetric matrix
        W=triu(W);      % ensure symmetry is preserved
        ud=2;           % halve number of removed links
    else
        ud=1;
    end
    
    ind=find(W);                            % find all links
    E=sortrows([ind W(ind)], -2);           % sort by magnitude
    en=round((Nroi^2-Nroi)*0.01*th/ud);     % number of links to be preserved. the value in front of '/ud' will be changed
    W(E(en+1:end,1))=0;                     % apply threshold
    
    % if symmetric, reconstruct symmetry
    if ud==2
        W=W+W.';
    end
    % binarize
    W=double(W~=0);
    binCONN(sidx, :, :) = W;
%     save([outpath, 'revision/1.binconn5/sub', pad(num2str(sidx, '%d'), 3, 'left', '0'), '.mat'], 'binconn');
end
save([outpath, 'revision/control_binCONN.mat'], 'binCONN', '-v7.3')

%% 4) SFC analysis
Nstep = 7;
SFC = zeros(Nsub, Nroi, Nroi, Nstep);
for sidx = 1 : Nsub
    disp(['subject = ', num2str(sidx)])
%     load([outpath, 'revision/1.binconn5/sub', pad(num2str(sidx, '%d'), 3, 'left', '0'), '.mat'])
    [sfc, ~, ~] = findwalks(squeeze(binCONN(sidx,:,:)));
    sfc = sfc(:, :, 1:Nstep);
    for step = 1 : Nstep
        W = reshape(sfc(:,:,step), [Nroi, Nroi]);
        sfc(:, :, step) = W - diag(diag(W));
    end
    sfc = normalize_sfc(sfc);
    SFC(sidx, :, :, :) = sfc;
%     save([outpath, 'revision/2.sfc5/sub', pad(num2str(sidx, '%d'), 3, 'left', '0'), '.mat'], 'sfc');
end
save([outpath, 'revision/control_norm_SFC.mat'], 'SFC', '-v7.3')

%% 5) Compute group average SFC
seed_idx = [7, 8, 9, 20, 27, 28, 54, 58, 59, 60, 89, 93, 99, 109, 152, 162, 165, 198, 212, 222, 231, 233, 234, 237, 238, 243]';

grpmean_SFC = cell(2,1);
grpmean_SFC{1} = zeros(length(seed_idx), Nroi, Nstep);
grpmean_SFC{2} = zeros(length(seed_idx), Nroi, Nstep);
for sidx = 1 : Nsub
    if group(sidx) == 1
%         load([inpath, '2.sfc5/sub', pad(num2str(sidx, '%d'), 3, 'left', '0'), '.mat']);
        grpmean_SFC{1} = grpmean_SFC{1} + squeeze(SFC(sidx,seed_idx,:,:));
    elseif group(sidx) == 2
%         load([inpath, '2.sfc5/sub', pad(num2str(sidx, '%d'), 3, 'left', '0'), '.mat']);
        grpmean_SFC{2} = grpmean_SFC{2} + squeeze(SFC(sidx,seed_idx,:,:));
    end
end
grpmean_SFC{1} = grpmean_SFC{1} / sum(group == 1);
grpmean_SFC{2} = grpmean_SFC{2} / sum(group == 2);

grpmean_DC = zeros(2, Nroi, Nstep);
grpmean_DC(1,:,:) = squeeze(sum(grpmean_SFC{1}, 1));
grpmean_DC(2,:,:) = squeeze(sum(grpmean_SFC{2}, 1));
save([outpath, 'revision/control_groupmeanSFC.mat'], 'grpmean_SFC', 'grpmean_DC')

%% 6) Calculate DC vector for all subjects
load([inpath, 'cluster_Fan_Net_r280.mat'])
net8 = cluster_Fan_Net.dat(1:246, 3);
net = cluster_Fan_Net.descrip{3,2};
net(9) = [];
num_network = 8;

seed_idx = [7, 8, 9, 20, 27, 28, 54, 58, 59, 60, 89, 93, 99, 109, 152, 162, 165, 198, 212, 222, 231, 233, 234, 237, 238, 243]';
load([outpath, 'a_group.mat'])
roi_dc = zeros(Nsub, Nroi, Nstep);
net_dc = zeros(Nsub, Nstep, num_network);
for sidx = 1 : Nsub
%     load([inpath, '2.sfc5/sub', pad(num2str(sidx, '%d'), 3, 'left', '0'), '.mat'])
    sfc = squeeze(SFC(sidx,seed_idx,:,:));
    for step = 1 : Nstep
        dc = sum(sfc(:,:,step), 1);
        dc(isinf(dc)|isnan(dc)) = 0;
        roi_dc(sidx, :, step) = dc;
        for nidx = 1 : 8
            net_dc(sidx, step, nidx) = mean(dc(net8 == nidx));
        end
    end
end
save([outpath, 'revision/control_wholesub_norm_roidc.mat'], 'roi_dc');
save([outpath, 'revision/control_wholesub_norm_netdc.mat'], 'net_dc');

%% 7) Plot hub regions
% Surface for 'plot_hemispheres'
parc = zeros(64984,1);
gii = gifti([inpath, 'conte69/BNA_conte69.L.func.gii']);
parc(1:32492,:) = gii.cdata;
gii = gifti([inpath, 'conte69/BNA_conte69.R.func.gii']);
parc(1+32492:end,:) = gii.cdata;
[surf_lh, surf_rh] = load_conte69();

group_name = {'HW', 'OB'};
load([inpath, 'cool_warm.mat'])
cool_warm(1,:) = cool_warm(129,:);

load([outpath, 'revision/control_groupmeanSFC.mat'])
clear grpmean_SFC
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
        saveas(gcf, [outpath, 'revision/figures/control/hub_', group_name{gidx}, num2str(mode), '.png'])
        close(gcf)
    end
end

% subcortex  --  only found in step 1
for gidx = 1 : 2
    for step = 1
        hub_1step = hub{gidx}(:, step);
        plot_subcortical(hub_1step(211:end), 'ventricles', 'False', 'cmap', 'RdBu_r', 'color_range', [-1.3 1.3])
        enigma_colormap(cool_warm)
        saveas(gcf, [outpath, 'revision/figures/control/hub_subcor_', group_name{gidx}, num2str(step), '.png'])
    end
end

%% 8) Group difference: roi-level
load([outpath, 'revision/control_wholesub_norm_roidc.mat']);
Nroi = 224;
H = zeros(Nroi, Nstep);
P = zeros(Nroi, Nstep);
T = zeros(Nroi, Nstep);

for step = 1 : Nstep
    ob_dc = [roi_dc(group==2,1:210,step), mean_subcor(roi_dc(group==2,:,step)')'];
    hw_dc = [roi_dc(group==1,1:210,step), mean_subcor(roi_dc(group==1,:,step)')'];
    for roi = 1 : Nroi
%         [~,p,~,stats] = ttest2(roi_dc(group==2,roi,step), roi_dc(group==1,roi,step));
        [~,p,~,stats] = ttest2(ob_dc(:,roi), hw_dc(:,roi));
        P(roi, step) = p;
        T(roi, step) = stats.tstat;
    end
    [selected, ~, ~, corrected_P] = fdr_bh(P(:, step), 0.05);
    H(:, step) = selected;
    P(:, step) = corrected_P;
end
save([outpath, 'revision/control_roi_ttest.mat'], 'H', 'P', 'T');

% Plot results
sigT = H .* T;
steps = [[1, 2, 3, 4]; [5, 6, 7, 7]];
clim = [[-5 -5 -5 -5]; [5 5 5 5]];
for mod = 1 : 2
    obj = plot_hemispheres2(sigT(1:210,steps(mod,:)), {surf_lh,surf_rh}, 'parcellation', parc, 'clim', clim);
    colormap(obj.figure, cool_warm)
    saveas(gcf, [outpath, 'revision/figures/control/ROIttest', num2str(mod), '.png'])
    close(gcf)
end

% subcortex
for gidx = 1 : 2
    for step = 1 : 7
        plot_subcortical(sigT(211:end, step), 'ventricles', 'False', 'cmap', 'RdBu_r', 'color_range', [-5 5])
        enigma_colormap(cool_warm)
        saveas(gcf, [outpath, 'revision/figures/control/ROIttest_subcor_', group_name{gidx}, num2str(step), '.png'])
%         close(gcf)
    end
end

%% 9) Group difference: network-level
load([outpath, 'revision/control_wholesub_norm_netdc.mat']);
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
    saveas(gcf, [outpath, 'revision/figures/control/NETttest_', num2str(nidx),'.png']);
end
%% 7) Plot DC
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
        saveas(gcf, [outpath, 'revision/figures/control/dc_', group_name{gidx}, num2str(mod), '.png'])
        close(gcf)
    end
end

% subcortex
for gidx = 1 : 2
    for step = 1 : Nstep
        dc_subcor = mean_subcor(squeeze(grpmean_DC(gidx, :, step))');
        plot_subcortical(dc_subcor, 'ventricles', 'False', 'cmap', 'RdBu_r', 'color_range', dc_range(step,gidx*2-1:gidx*2));%clim{gidx}(:,step)')
        enigma_colormap(cool_warm)
        saveas(gcf, [outpath, 'revision/figures/control/dc_subcor_', group_name{gidx}, num2str(step), '.png'])
%         close(gcf)
    end
end
end