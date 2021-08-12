function perm_clinical()
%% R.11) Bootstrapping
% (1) Report probability map for seed regions
% (2) Report probability map for SFC DC values
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

%% Permutaion - R
% eat_score = xlsread([store7, 'hblee/2.stepwise/Enhanced_NKI.xlsx'], 9, 'H3:N303'); %[edeq-4; tfeq-3]
% save([outpath, 'a_eating_scores.mat'], 'eat_score');
load([outpath, 'a_eating_scores.mat']);
tfeq = eat_score(:, 5:end);
total_tfeq = sum(tfeq, 2);
tfeq = [tfeq, sum(tfeq(:,1:2), 2), sum(tfeq, 2)];
clear eat_score total_tfeq
%%
load([outpath, 'd_wholesub_norm_netdc.mat']);
% network = [1, 2, 3, 5, 8];
network = [1, 2, 3, 4, 5, 6, 7, 8];
netdc = squeeze(net_dc(1,:,7,network));
Nperm = 5000;
num_score = size(tfeq, 2);
num_network = length(network);
perm_R = zeros(Nperm, num_network, num_score);
% perm_P = zeros(Nperm, num_network, num_score);

perm = 0;
while perm < Nperm
%     if rem(perm+1, 500) == 0
%         disp(perm+1);
%     end
    tfeq_perm = tfeq(randperm(Nsub), :);
    for score = 1 : num_score
        for nidx = 1 : num_network
            [r, ~] = corrcoef(tfeq_perm(:,score), netdc(:,nidx));
            perm_R(perm+1, nidx, score) = r(1,2);
%             perm_P(perm+1, nidx, score) = p(1,2);
        end
    end
    perm = perm + 1;
end
save([outpath, 'revision/tfeq_perm_corr.mat'], 'perm_R')

%% real R
real_R = zeros(num_network, num_score);
th = zeros(num_network, num_score, 2);
for score = 1 : num_score
    for nidx = 1 : num_network
        [r, ~] = corrcoef(tfeq(:,score), netdc(:,nidx));
        real_R(nidx, score) = r(1,2);
        Rs = sort(perm_R(:,nidx, score));
        th(nidx, score, 1) = Rs(Nperm*0.025);
        th(nidx, score, 2) = Rs(Nperm - Nperm*0.025 + 1);
%         histogram(perm_R(:,1,1), 'facecolor',[0.5 0.5 0.5])
%         saveas(gcf, [outpath, 'figures/norm/perm/tfeq', num2str(score), '-net', num2str(nidx), '.png'])
    end
end

%% plot histogram
P = zeros(num_network, num_score);
for score = 1 : num_score
    for nidx = 1 : num_network
        [r, ~] = corrcoef(tfeq(:,score), netdc(:,nidx));
        real_R = r(1,2);
        if real_R > 0
            P(nidx, score) = length(find(perm_R(:,nidx, score) > real_R)) / Nperm * 2;
        else
            P(nidx, score) = length(find(perm_R(:,nidx, score) < real_R)) / Nperm * 2;
        end
%         histogram(perm_R(:,1,1), 'facecolor',[0.5 0.5 0.5])
%         saveas(gcf, [outpath, 'figures/norm/perm/tfeq', num2str(score), '-net', num2str(nidx), '.png'])
    end
end

%% FDR correction
% 1) across all tests
% [H, ~, ~, P_corrected] = fdr_bh(P, 0.05);
% 2) across all tests, only for networks with
%    significant group difference in step 7
network = [1, 2, 3, 5, 8];
[H, ~, ~, P_corrected] = fdr_bh(P(network,:), 0.05);
% 3) across scores
% P_corrected = zeros(num_network, num_score);
% for nidx = 1 : num_network
%     [~, ~, ~, p] = fdr_bh(P(nidx,:), 0.05);
%     P_corrected(nidx, :) = p;
% end
% H = P_corrected < 0.05;
% 4) across networks, only for networks with
%    significant group difference in step 7
% network = [1, 2, 3, 5, 8];
% P_corrected = zeros(5, num_score);
% for score = 1 : num_score
%     [~, ~, ~, p] = fdr_bh(P(network, score), 0.05);
%     P_corrected(:, score) = p;
% end
% H = P_corrected < 0.05;

%%
R = zeros(num_network, num_score);
for score = 1 : num_score
    for nidx = 1 : num_network
        [r, ~] = corrcoef(tfeq(:,score), netdc(:,nidx));
        R(nidx, score) = r(1,2);
    end
end

end