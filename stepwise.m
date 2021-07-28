%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Perform stepwise connectivity analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stepwise()
if ispc
    store7 = 'V:/';
else
    store7 = '/store7/';
end

inpath = [store7, 'hblee/2.stepwise/8.regress-466/'];
outpath = [store7, 'hblee/2.stepwise/9.regress-301/'];
Nsub = 301;
Nroi = 246;
Nstep = 200;

%% Subject selection
issub = xlsread([store7, 'hblee/2.stepwise/Enhanced_NKI.xlsx'], 8, 'L2:L467')==1;
save([outpath, 'subj_selection.mat'])

load([inpath, 'a_dataset.mat'])
age = age(issub); 
sex = sex(issub);
obesity = obesity(issub, :);
group = group(issub);
subID = subID{issub};
save([outpath, 'a_dataset.mat'], 'age', 'sex', 'obesity', 'group')

group = xlsread([store7, 'hblee/2.stepwise/Enhanced_NKI.xlsx'], 9, 'G2:G302');
save([outpath, 'a_bingroup.mat'], 'group')

load([inpath, 'c_headmotion.mat'])
headmotion = headmotion(issub);
save([outpath, 'a_headmotion.mat'], 'headmotion')

load([inpath, '1.ridge/b_CONN_raw.mat'])
CONN = CONN(issub, : , :);
save([outpath, 'a_conn_ridge.mat'], 'CONN')

%% Regress out from BMI & WHR
% make model
tbl = table(obesity(:,1), age, sex, 'VariableNames', {'bmi','age','sex'});
mdl = fitlm(tbl, 'bmi ~ age+sex');
bmi_reg = mdl.Residuals.Raw;
tbl = table(obesity(:,2), age, sex, 'VariableNames', {'whr','age','sex'});
mdl = fitlm(tbl, 'whr ~ age+sex');
whr_reg = mdl.Residuals.Raw;
save([outpath, 'a_obesity_reg.mat'], 'bmi_reg', 'whr_reg')

%% Save connectivity matrix
for v = 1 : 3
    CONN = zeros(Nsub, Nroi, Nroi);
    for i = 1 : Nsub
        load([basepath{v}, '1.conn/sub', pad(num2str(i, '%d'), 3, 'left', '0'), '.mat']);
        CONN(i,:,:) = conn;
    end
    save([basepath{v}, 'b_CONN_raw.mat'], 'CONN', '-v7.3');
end

%% Regress out age & sex & headmotion from CONNECTIVITY
% load([outpath, 'a_dataset.mat'])
% load([outpath, 'a_headmotion.mat'])
load([outpath, 'a_conn_ridge.mat'])
CONN_reg = zeros(Nsub, Nroi, Nroi);
for i = 1 : Nroi
    for j = 1 : Nroi
        if i ~= j
            % make model
            tbl = table(CONN(:,i,j), age, sex, headmotion, 'VariableNames', {'R','age','sex', 'headmotion'});
            mdl = fitlm(tbl, 'R ~ age+sex+headmotion');
            CONN_reg(:,i,j) = mdl.Residuals.Raw;
        end
    end
end
save([outpath, 'b_conn_regressed2.mat'], 'CONN_reg', '-v7.3')

%% Save connectivity image
for i = 1 : Nsub
    conn = squeeze(CONN_reg(i,:,:));
    conn_lim = max(max(conn(:)), -min(conn(:)));
    imagesc(conn), colorbar, axis off, caxis([-conn_lim conn_lim])
    saveas(gcf, [outpath, '1.conn2/sub', pad(num2str(i, '%d'), 3, 'left', '0'), '.png']);
end
%%
for i = 1 : Nsub
    for v = 1 : 3
        subplot(3,2,2*v-1), imagesc(squeeze(raw{v}(i,:,:))), colorbar, axis off
        subplot(3,2,2*v), imagesc(squeeze(regressed{v}(i,:,:))), colorbar, axis off
        saveas(gcf, [outpath, 'figure/conn/sub', pad(num2str(i, '%d'), 3, 'left', '0'), '.png']);
    end
end
%% Find WHR-correlated ROIs
v = 1;
load([outpath, 'a_dataset.mat'])
% load([outpath, 'b_obesity_reg.mat'])
load([outpath, 'a_conn_ridge.mat'])
DC = sum(CONN, 3);
for roi = 1 : Nroi
    [r, p] = corrcoef(obesity(:,1), DC(:,roi));
%     [r, p] = corrcoef(bmi_reg, DC(:,roi));
    R(roi, 1) = r(1,2);
    P(roi, 1) = p(1,2);
end
[selected, ~, ~, corrected_P] = fdr_bh(P, 0.05);
seed_idx = find(selected == 1);
significant_R = R(seed_idx);
significant_P = corrected_P(seed_idx);
buf = corrected_P .* (corrected_P<=0.05);
disp(['p1 = ', num2str(sum(P<0.05)), ' & p2 = ', num2str(length(seed_idx))])

%% Visualize seed on surface
bna = load_nii('Z:/hblee/ETC/Atlas/BNA/BNA_3mm.nii');
roi = bna.img;

WHR_related = zeros(size(roi));
for i = 1 : length(seed_idx)
    idx = find(roi(:) == seed_idx(i));
    WHR_related(idx) = significant_R(i);
end

bna.img = WHR_related;
save_nii(bna, [outpath, 'image/seed.nii']);

%% Binarize correlation matrix
disp(['Binarization of connectivity - processing', newline]);
opt_sparsity = zeros(Nsub, 1);
for sidx = 1 : Nsub
    disp(['subject = ', num2str(sidx)])
    conn = squeeze(CONN_reg(sidx, :, :));
    ConnMatZ_thr = zeros(Nroi, Nroi, 100);
    for th = 1:100
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
        ConnMatZ_thr(:,:,th) = W(:,:);
    end
    
    % Check the optimal sparsity
    for th = 1:100
        adj = ConnMatZ_thr(:,:,th);
        % Stop when all indices becomes 1
        num = get_components(adj);
        if min(num) == 1 && max(num) == 1
            break
        end
    end
    opt_sparsity(sidx) = th;  % the minimum percentage(th) so that the network is fully coinnected
    binconn = ConnMatZ_thr(:,:,opt_sparsity(sidx));
    save([outpath, '2.binconn/sub', pad(num2str(sidx, '%d'), 3, 'left', '0'), '.mat'], 'binconn');
end
save([outpath, 'b_binTH.mat'], 'opt_sparsity')
figure, histogram(opt_sparsity), title('distribution of binarization threshold')
saveas(gcf, [outpath, 'figure/hist-binTH.png'])

%% Fix binarization threshold
disp(['Binarization of connectivity - processing', newline]);
opt_sparsity = zeros(Nsub, 1);
th = 5;
for sidx = 1 : Nsub
    disp(['subject = ', num2str(sidx)])
    % load([basepath, '1.conn0.5/conn-sub', pad(num2str(sidx, '%d'), 3, 'left', '0'), '.mat'])
    conn = squeeze(CONN(sidx, :, :));
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
    binconn = W;
    save([outpath, '1.binconn5/sub', pad(num2str(sidx, '%d'), 3, 'left', '0'), '.mat'], 'binconn');
end

%% SFC analysis
Nstep = 200;
for sidx = 1 : Nsub
    disp(['subject = ', num2str(sidx)])
    load([outpath, '1.binconn5/sub', pad(num2str(sidx, '%d'), 3, 'left', '0'), '.mat'])
    [sfc, ~, ~] = findwalks(binconn);
    sfc = sfc(:, :, 1:Nstep);
    for step = 1 : Nstep
        W = reshape(sfc(:,:,step), [Nroi, Nroi]);
        sfc(:, :, step) = W - diag(diag(W));
    end
    save([outpath, '2.sfc5/sub', pad(num2str(sidx, '%d'), 3, 'left', '0'), '.mat'], 'sfc');
end

%% Compute group average SFC
seed_idx = [7, 8, 9, 20, 27, 28, 54, 58, 59, 60, 89, 93, 99, 109, 152, 162, 165, 198, 212, 222, 231, 233, 234, 237, 238, 243]';
seed_sub{1} = [7, 8, 20, 27, 28, 89, 93, 99, 109, 152, 198]';   % Cognition
seed_sub{2} = [165, 237, 238, 243]';                            % Reward
seed_sub{3} = [9, 54, 58, 59, 60, 162, 222, 231, 233, 234]';    % Sensorimotor

seed_avg = cell(2,1);
seed_avg{1} = zeros(length(seed_idx), Nroi, Nstep);
seed_avg{2} = zeros(length(seed_idx), Nroi, Nstep);
for sidx = 1 : Nsub
    load([outpath, '2.sfc5/sub', pad(num2str(sidx, '%d'), 3, 'left', '0'), '.mat']);
    if group(sidx) == 1
        seed_avg{1} = seed_avg{1} + sfc(seed_idx,:,:);
    elseif group(sidx) == 2
        seed_avg{2} = seed_avg{2} + sfc(seed_idx,:,:);
    end
end
seed_avg{1} = seed_avg{1} / sum(group == 1);
seed_avg{2} = seed_avg{2} / sum(group == 2);
save([outpath, 'c_seed_avg.mat'], 'seed_avg')

%% Save seed connectivity figure matrix
group_name = {'HW', 'OB'};
for gidx = 1 : 2
    for step = [1, 2, 3, 5, 7, 8, 10, 20, 30, 40, 50, 60, 80, 100, 150, 200]
        imagesc(seed_avg{gidx}(:,:,step)), axis off %, colorbar
        yticks([1:length(seed_idx)]), yticklabels(seed_idx), xticklabels([])
        saveas(gcf, [outpath, 'figure/sfc/', group_name{gidx}, '-step', num2str(step), '.png']);
    end
end

%% Save DC image
bna_path = 'Z:/hblee/ETC/Atlas/BNA/BNA_3mm.nii';
bna = load_nii(bna_path);
roi = bna.img;
% For rich/local nodes
group_name = {'HW', 'OB'};

load([outpath, 'c_seed_avg.mat'])
for gidx = 1 : 2
    seedconn = seed_avg{gidx};
    for step = [2] %[1, 3, 5, 10, 25, 50, 75, 100, 200]
        image = zeros(size(roi), 'like', roi);
        dc = rescale(sum(seedconn(:,:,step), 1));
        dc(isinf(dc)|isnan(dc)) = 0;
        for i = 1 : Nroi
            idx = find(roi(:) == i);
            image(idx) = dc(i);
        end
        bna.img = image;
        save_nii(bna, [outpath, 'image/dc/', group_name{gidx}, '-step', num2str(step),'.nii']);
    end
end

%% Save dc figure
surf_file = 'Z:\hblee\Toolbox\BrainNetViewer\Data\SurfTemplate\BrainMesh_ICBM152_smoothed.nv';
opt_file = 'Z:\hblee\2.stepwise\6.refinement\option_dc.mat';
loadpath = [outpath, 'image\dc\'];
savepath = [outpath, 'figure\dc\'];
for gidx = 1 : 2
    for step = [2] %[1, 3, 5, 10, 25, 50, 75, 100, 200]
        if gidx == 1
            disp([newline,'Healthy - step',num2str(step),newline])
            vol_file = [loadpath, 'HW-step',num2str(step),'.nii'];
            fig_file = [savepath, 'HW-',num2str(step),'.jpg'];
        else
            disp([newline,'Obesity - step',num2str(step),newline])
            vol_file = [loadpath, 'OB-step',num2str(step),'.nii'];
            fig_file = [savepath, 'OB-',num2str(step),'.jpg'];
        end
        BrainNet_MapCfg(surf_file, vol_file, opt_file, fig_file);
    end
end
disp([newline,'Finish saving images', newline])

%% Spider plot of R and P
load([outpath, 'cluster_Fan_Net_r280.mat'])
net8 = cluster_Fan_Net.dat(1:246, 3);
net = cluster_Fan_Net.descrip{3,2};
net(9) = [];
seed_net = net8(seed_idx);
for nidx = 1 : 8
%     net_R(nidx) = mean(significant_R(seed_net == nidx));
%     net_P(nidx) = mean(significant_P(seed_net == nidx));
    net_R(nidx) = mean(R(net8 == nidx));
    net_P(nidx) = mean(corrected_P(net8 == nidx));
end
% spider_plot([net_R; net_P], 'AxesLabels', net);
% legend('r-value', 'p-value', 'Location', 'southoutside');
disp(['min: ', num2str(min(net_R)), ' - max: ', num2str(max(net_R))])
disp(['min: ', num2str(min(net_P)), ' - max: ', num2str(max(net_P))])
% spider_plot2(net_R, 'AxesLabels', net, 'AxesLimits', [repmat(-0.2, 1, 8); repmat(0.2, 1, 8)], 'AxesInterval', 4);
spider_plot2(net_R, 'AxesLabels', net, 'AxesLimits', [repmat(-0.1, 1, 8); repmat(0.1, 1, 8)], 'AxesInterval', 4);
saveas(gcf, [outpath, 'figure/r-value_whole.png']);
% spider_plot2(net_P, 'AxesLabels', net, 'AxesLimits', [repmat(0.0, 1, 8); repmat(0.05, 1, 8)], 'AxesInterval', 5);
spider_plot2(net_P, 'AxesLabels', net, 'AxesLimits', [repmat(0.0, 1, 8); repmat(0.6, 1, 8)], 'AxesInterval', 4);
saveas(gcf, [outpath, 'figure/p-value_whole.png']);

%% Spider plot of DC
load([outpath, 'cluster_Fan_Net_r280.mat'])
net8 = cluster_Fan_Net.dat(1:246, 3);
net = cluster_Fan_Net.descrip{3,2};
net(9) = [];
seed_net = net8(seed_idx);
a = zeros(2, Nroi);

load([outpath, 'c_seed_avg.mat'])
for step = 2 % [1, 3, 5, 10, 20, 25, 40, 50, 60, 75, 80, 100, 150, 200]
    for gidx = 1 : 2
        seedconn = seed_avg{gidx};
        dc = sum(seedconn(:,:,step), 1);
        dc(isinf(dc)|isnan(dc)) = 0;
        % a(gidx,:) = dc;
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
    saveas(gcf, [outpath, 'figure/spider/', num2str(v), '-step', num2str(step),'.png']);
end
end
