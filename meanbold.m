%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Average BOLD signal for each atlas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function meanbold()
basepath = 'X:/path/myfolder/';
%% 1) Save demogrphic information
sex = xlsread([basepath, 'Enhanced_NKI.xlsx'], 6, 'C2:C469');
obesity = xlsread([basepath, 'Enhanced_NKI.xlsx'], 6, 'D2:E469');
group = xlsread([basepath, 'Enhanced_NKI.xlsx'], 6, 'F2:F469');
save([basepath, 'a_dataset.mat'], 'sex', 'obesity', 'group')

load([basepath, 'a_sublist.mat'])
Nsub = length(sublist);
filepath = 'REST_645/func_results_REST_645/Smooth_REST_645.nii.gz';

%% 2) Get BNA atlas map
bna_file = [basepath, 'Atlas/BNA_3mm.nii'];
bna = load_nii(bna_file);
atlas = bna.img;
Nroi = max(atlas(:));

%% 3) Average bold signal for each ROI
for sidx = 1 : Nsub
    subID = sublist{sidx};
    disp(strcat(['Subject ',int2str(sidx),' -- ', subID]));
    datapath = [basepath, 'data/'];
        disp('ERROR: subject folder not exists');
    end
    % Load fMRI data
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
end
