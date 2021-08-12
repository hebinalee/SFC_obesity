%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Compute stepwise connectivity matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stepwise_fc()
inpath = 'X:/path/myfolder/inputs/';
outpath = 'X:/path/myfolder/outputs/';
Nsub = 301;
Nroi = 246;
Nstep = 200;

%% 1) Binarize connectivity matrix with threshold
disp(['## Binarization of connectivity - processing', newline]);
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


%% 2) Construct SFC matrix
Nstep = 200;
for sidx = 1 : Nsub
    disp(['subject = ', num2str(sidx)])
    load([outpath, '1.binconn5/sub', pad(num2str(sidx, '%d'), 3, 'left', '0'), '.mat'])
    [sfc, ~, ~] = findwalks(binconn);
    sfc = sfc(:, :, 1:Nstep);
    for step = 1 : Nstep
        W = reshape(sfc(:,:,step), [Nroi, Nroi]);
        W = W - diag(diag(W));                           % off-diagonal
        sfc(:, :, step) = (W - mean(W(:))) / std(W(:));  % Normalization
    end
    save([outpath, '2.sfc5/sub', pad(num2str(sidx, '%d'), 3, 'left', '0'), '.mat'], 'sfc');
end


%% 3) Save DC values per ROI/Network for all subjects
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
    load([inpath, '2.sfc5/sub', pad(num2str(sidx, '%d'), 3, 'left', '0'), '.mat'])
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
end
