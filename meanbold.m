%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Average BOLD signal for each atlas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function meanbold(a)
basepath = 'V:/hblee/2.stepwise/6.refinement/';
sex = xlsread('V:/hblee/2.stepwise/Enhanced_NKI.xlsx', 6, 'C2:C469');
obesity = xlsread('V:/hblee/2.stepwise/Enhanced_NKI.xlsx', 6, 'D2:E469');
group = xlsread('V:/hblee/2.stepwise/Enhanced_NKI.xlsx', 6, 'F2:F469');
save([basepath, 'a_dataset.mat'], 'sex', 'obesity', 'group')
% 
% histogram(obesity(:,1),[2.5:2.5:50]), title('distribution of BMI')
% hold on, line([25 25],ylim,'LineStyle','--')
% hold on, line([30 30],ylim,'LineStyle','--')
% histogram(obesity(sex==1,2),[0.55:0.05:1.3]), title('distribution of WHR for male'), hold on, line([0.9 0.9],ylim,'LineStyle','--')
% histogram(obesity(sex==2,2),[0.55:0.05:1.3]), title('distribution of WHR for female'), hold on, line([0.85 0.85],ylim,'LineStyle','--')

load([basepath, 'a_sublist.mat'])
Nsub = length(sublist);
filepath = 'REST_645/func_results_REST_645/Smooth_REST_645.nii.gz';

% Get BNA atlas map
bna_file = 'V:/hblee/ETC/Atlas/BNA/BNA_3mm.nii';
bna = load_nii(bna_file);
atlas = bna.img;
Nroi = max(atlas(:));

for sidx = a : Nsub
    subID = sublist{sidx};
    disp(strcat(['list = ',int2str(sidx),' -- ', subID]));
    % calculate mean BOLD signal of region
    datapath1 = 'U:/hblee/1.sCCA/eNKI/';
    datapath2 = 'V:/hblee/2.stepwise/eNKI-process/eNKI_add/';
    if isdir([datapath1, subID])
        datapath = datapath1;
        is(sidx) = 1;
    elseif isdir([datapath2, subID])
        datapath = datapath2;
        is(sidx) = 1;
    else
        disp('ERROR: subject folder not exists');
    end
    %% (1) Load fMRI data
    sbj_file = fullfile(datapath, subID, filepath);
    sbj = load_nii(sbj_file);
    fmri = sbj.img;
    if size(fmri, 1) == 61
        fmri = fmri(2:61,2:73,2:61, :);
    end
    % Prepare time series and calculate mean BOLD signal
    fmri = reshape(fmri, [length(atlas(:)), size(fmri,4)]);
    for j = 1 : Nroi
        region_idx = find(atlas(:) == j);     % find voxels included in the region
        region_bold = fmri(region_idx, :);    % find BOLD signal of those voxels
        bold(:, j) = mean(region_bold, 1);     % calculate mean BOLD signal of region (time point x NRoi)
    end
    save([basepath, '0.meanBOLD/bold-sub', pad(num2str(sidx, '%d'), 3, 'left', '0'), '.mat'], 'bold');
    clear bold
end
