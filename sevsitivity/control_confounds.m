function control_confounds()
inpath = 'X:/path/myfolder/inputs/';
outpath = 'X:/path/myfolder/outputs/';
Nsub = 301;
Nroi = 246;
Nstep = 7;


%% 1) Regress out from BMI & WHR
make model
tbl = table(obesity(:, 1), age, sex, 'VariableNames', {'bmi','age','sex'});
mdl = fitlm(tbl, 'bmi ~ age+sex');
bmi_reg = mdl.Residuals.Raw;
tbl = table(obesity(:, 2), age, sex, 'VariableNames', {'whr','age','sex'});
mdl = fitlm(tbl, 'whr ~ age+sex');
whr_reg = mdl.Residuals.Raw;
save([outpath, 'obesity_regressed.mat'], 'bmi_reg', 'whr_reg')

%% 2) Regress out age & sex & headmotion from CONNECTIVITY
% load([outpath, 'a_dataset.mat'])
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
save([outpath, 'CONN_regressed.mat'], 'CONN_reg', '-v7.3')
end
