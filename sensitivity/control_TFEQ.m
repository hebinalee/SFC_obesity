function control_confounds()
inpath = 'X:/path/myfolder/inputs/';
outpath = 'X:/path/myfolder/outputs/';
Nsub = 301;
Nroi = 246;
Nstep = 7;


%% 1) Regress out TFEQ scores from BMI & WHR
disp(['## Regress out TFEQ from obesity phenotypes - processing', newline]);
load([outpath, 'a_dataset.mat'])
load([outpath, 'a_eating_scores.mat'])

tbl = table(obesity(:, 1), tfeq(:, 1), tfeq(:, 2), 'VariableNames', {'bmi','tfeq1','tfeq2'});
mdl = fitlm(tbl, 'bmi ~ tfeq1 + tfeq2');
bmi_reg = mdl.Residuals.Raw;
tbl = table(obesity(:, 2), tfeq(:, 1), tfeq(:, 2), 'VariableNames', {'whr','tfeq1','tfeq2'});
mdl = fitlm(tbl, 'whr ~ tfeq1 + tfeq2');
whr_reg = mdl.Residuals.Raw;
save([outpath, 'controlTFEQ/obesity_regressed_out.mat'], 'bmi_reg', 'whr_reg')


%% 2) Regress out TFEQ scores from CONNECTIVITY
disp(['## Regress out TFEQ from connectivity - processing', newline]);
load([inpath, 'a_conn_ridge.mat'])
CONN_reg = zeros(Nsub, Nroi, Nroi);
for i = 1 : Nroi
    disp(['i = ', num2str(i)])
    for j = 1 : Nroi
        if i ~= j
            % make model
            tbl = table(CONN(:,i,j), tfeq(:, 1), tfeq(:, 2), 'VariableNames', {'R','tfeq1','tfeq2'});
            mdl = fitlm(tbl, 'R ~ tfeq1 + tfeq2');
            CONN_reg(:,i,j) = mdl.Residuals.Raw;
        end
    end
end
save([outpath, 'controlTFEQ/CONN_regressed_out.mat'], 'CONN_reg', '-v7.3')
end


%% 3) Find WHR-correlated ROIs (DEFINE SEED)
disp(['## Define seed regions - processing', newline]);
DC = sum(CONN_reg, 3);
for roi = 1 : Nroi
    [r, p] = corrcoef(whr_reg, DC(:,roi));
    R(roi, 1) = r(1,2);
    P(roi, 1) = p(1,2);
end
[selected, ~, ~, corrected_P] = fdr_bh(P, 0.05);
seed_idx = find(selected == 1);
significant_R = R(seed_idx);
significant_P = corrected_P(seed_idx);
R_seed = selected .* R;
disp(['p1 = ', num2str(sum(P<0.05)), ' & p2 = ', num2str(length(seed_idx))])
save([outpath, 'controlTFEQ/Rvalue_WHRrelated.mat'], 'R', 'R_seed')


%% 4) Binarize connectivity matrix with threshold 5
disp(['## Binarization of connectivity - processing', newline]);
th = 5;
for sidx = 1 : Nsub
    disp(['subject = ', num2str(sidx)])
    conn = squeeze(CONN(sidx, :, :));
    binconn = binarize_conn(conn);
    save([outpath, 'controlTFEQ/1.binconn/sub', pad(num2str(sidx, '%d'), 3, 'left', '0'), '.mat'], 'binconn');
end


%% 5) Construct SFC matrix
Nstep = 5;
for sidx = 1 : Nsub
    disp(['subject = ', num2str(sidx)])
    load([outpath, 'controlTFEQ/1.binconn/sub', pad(num2str(sidx, '%d'), 3, 'left', '0'), '.mat'])
    sfc = compute_sfc(binconn, Nstep)
    save([outpath, 'controlTFEQ/2.sfc/sub', pad(num2str(sidx, '%d'), 3, 'left', '0'), '.mat'], 'sfc');
end


%% 7) Save DC values per ROI/Network for all subjects
load([inpath, 'cluster_Fan_Net_r280.mat'])
net8 = cluster_Fan_Net.dat(1:246, 3);
net = cluster_Fan_Net.descrip{3,2};
net(9) = [];
num_network = 8;

load([inpath, 'c_seed_regions.mat']);
load([outpath, 'a_group.mat'])

roi_dc = zeros(Nsub, Nroi, Nstep);
net_dc = zeros(Nsub, Nstep, num_network);
for sidx = 1 : Nsub
    load([inpath, '2.sfc/sub', pad(num2str(sidx, '%d'), 3, 'left', '0'), '.mat'])
    for step = 1 : Nstep
        dc = sum(sfc(:,:,step), 1);
        dc(isinf(dc)|isnan(dc)) = 0;
        roi_dc(sidx, :, step) = dc;
        for nidx = 1 : 8
            net_dc(sidx, step, nidx) = mean(dc(net8 == nidx));
        end
    end
end
save([outpath, 'wholesub_ROI_dc.mat'], 'roi_dc');
save([outpath, 'wholesub_NET_dc.mat'], 'net_dc');
load([outpath, 'a_group.mat'])
