%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Compute connectivity matrix from meanBOLD signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function connectivity(a)
basepath = 'V:/hblee/2.stepwise/6.refinement/';
load([basepath, 'a_sublist.mat'])
Nsub = length(sublist);

for sidx = a : Nsub
    subID = sublist{sidx};
    disp(strcat(['list = ',int2str(sidx),' -- ', subID]));
    load([basepath, '0.meanBOLD/bold-sub', pad(num2str(sidx, '%d'), 3, 'left', '0'), '.mat'])
    % calculate mean BOLD signal of region
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
end
