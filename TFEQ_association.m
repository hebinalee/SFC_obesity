%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Associate degree of SFC with TFEQ scores
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TFEQ_association()
inpath = 'X:/path/myfolder/inputs/';
outpath = 'X:/path/myfolder/outputs/';
Nsub = 301;
Nroi = 246;
Nstep = 7;

%% 1) Load TFEQ scores
% eat_score = xlsread([store7, 'hblee/2.stepwise/Enhanced_NKI.xlsx'], 9, 'H3:N303'); %[edeq-4; tfeq-3]
load([outpath, 'a_eating_scores.mat']);
tfeq = eat_score(:, 5:end);
total_tfeq = sum(tfeq, 2);
tfeq = [tfeq, sum(tfeq(:,1:2), 2), sum(tfeq, 2)];
clear eat_score total_tfeq


%% 2) Load degree of SFC for all subjects
load([outpath, 'wholesub_NET_dc.mat']);
network = [1, 2, 3, 4, 5, 6, 7, 8];
netdc = squeeze(net_dc(1,:,7,network));
Nperm = 5000;
num_score = size(tfeq, 2);
num_network = length(network);
perm_R = zeros(Nperm, num_network, num_score);


%% 3) Perform null distribution using randomized data
perm = 0;
while perm < Nperm
     if rem(perm+1, 500) == 0
         disp(perm+1);
     end
    tfeq_perm = tfeq(randperm(Nsub), :);
    for score = 1 : num_score
        for nidx = 1 : num_network
            [r, ~] = corrcoef(tfeq_perm(:,score), netdc(:,nidx));
            perm_R(perm+1, nidx, score) = r(1,2);
        end
    end
    perm = perm + 1;
end
save([outpath, 'tfeq_perm_corr.mat'], 'perm_R')


%% 4) Compute real r-values with real data
real_R = zeros(num_network, num_score);
th = zeros(num_network, num_score, 2);
for score = 1 : num_score
    for nidx = 1 : num_network
        [r, ~] = corrcoef(tfeq(:,score), netdc(:,nidx));
        real_R(nidx, score) = r(1,2);
        Rs = sort(perm_R(:,nidx, score));
        th(nidx, score, 1) = Rs(Nperm*0.025);
        th(nidx, score, 2) = Rs(Nperm - Nperm*0.025 + 1);
    end
end


%% 5) Compute p-values
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
    end
end
[H, ~, ~, P_corrected] = fdr_bh(P, 0.05);
save([outpath, 'tfeq_correlation_analysis.mat'], 'H', 'P', 'real_R')
end
