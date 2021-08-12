%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Compute connectivity matrix from meanBOLD signals
%%  and define seed regions (WHR-related regions)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function connectivity_seed()
basepath = 'X:/path/myfolder/';
load([basepath, 'a_sublist.mat'])
Nsub = length(sublist);
Nroi = 246;

%% 1) Construct functional connectivity matrix
CONN = zeros(Nsub, Nroi, Nroi);
for sidx = 1 : Nsub
    subID = sublist{sidx};
    disp(strcat(['list = ',int2str(sidx),' -- ', subID]));
    load([basepath, '0.meanBOLD/bold-sub', pad(num2str(sidx, '%d'), 3, 'left', '0'), '.mat'])
    % Compute using ridge regression
    bold = bold';
    Nroi = size(bold, 2);
	rho = 0.5;  % UKBiobank style
    grot = cov(bold);
    grot = grot/sqrt(mean(diag(grot).^2));
    grot = -inv(grot+rho*eye(Nroi));
    grot = (grot ./ repmat(sqrt(abs(diag(grot))),1,Nroi)) ./ repmat(sqrt(abs(diag(grot)))',Nroi,1);
    grot(eye(Nroi)>0)=0;
    netmats2 = reshape(grot,1,Nroi*Nroi);
    R_part = reshape(netmats2,Nroi,Nroi);
    conn = .5*log((1+R_part)./(1-R_part));
    save([basepath, '1.conn0.5/conn-sub', pad(num2str(sidx, '%d'), 3, 'left', '0'), '.mat'], 'conn');
    CONN(sidx, :, :) = conn;
end
save([basepath, 'a_conn_ridge.mat'], 'CONN')

%% 2) Define seed regions by associating degree with WHR
load([basepath, 'a_dataset.mat'])
DC = sum(CONN, 3);
for roi = 1 : Nroi
    [r, p] = corrcoef(obesity(:,1), DC(:,roi));
    R(roi, 1) = r(1,2);
    P(roi, 1) = p(1,2);
end
[selected, ~, ~, corrected_P] = fdr_bh(P, 0.05);
seed_idx = find(selected == 1);
significant_R = R(seed_idx);
significant_P = corrected_P(seed_idx);
buf = corrected_P .* (corrected_P<=0.05);
disp([num2str(length(seed_idx)), 'regions are selected as seed'])
save([basepath, 'c_seed_regions.mat']);
end
