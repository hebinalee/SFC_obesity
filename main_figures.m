%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Plot the results of group analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function main_figures()
inpath = 'X:/path/myfolder/inputs/';
outpath = 'X:/path/myfolder/outputs/';
Nsub = 301;
Nroi = 246;
Nstep = 7;

%% 0) Get ready for plot
parc = zeros(64984,1);
gii = gifti([inpath, 'conte69/BNA_conte69.L.func.gii']);
parc(1:32492,:) = gii.cdata;
gii = gifti([inpath, 'conte69/BNA_conte69.R.func.gii']);
parc(1+32492:end,:) = gii.cdata;
[surf_lh, surf_rh] = load_conte69();

load([inpath, 'cool_warm.mat'])
cool_warm(1,:) = cool_warm(129,:);

group_name = {'HW', 'OB'};


%% 1) R-values
load([basepath, 'c_Rvalue.mat']);

% cortex
obj = plot_hemispheres2(R(1:210,:), {surf_lh,surf_rh}, 'parcellation', parc, 'clim', [-0.3; 0.3]);
colormap(obj.figure, cool_warm)
saveas(gcf, [outpath, 'figures/seed_R.png'])
close(gcf)

% subcortical
plot_subcortical(mean_subcortical(R), 'ventricles', 'False', 'cmap', 'RdBu_r', 'color_range', [-0.3 0.3]);
enigma_colormap(cool_warm)
saveas(gcf, [outpath, 'figures/seed_R_subcor.png'])
close(gcf)


%% 2) Degree centrality
load([outpath, 'groupmeanSFC.mat'])
clear grpmean_SFC

% cortex
for i = 1 : Nstep
    dc_range(i,:)=[min(grpmean_DC(1,1:210,i)), max(grpmean_DC(1,1:210,i)), min(grpmean_DC(2,1:210,i)), max(grpmean_DC(2,1:210,i))];
end

steps = [[1, 2, 3, 4]; [5, 6, 7, 7]];
for mod = 1 : 2         % low/high steps
    for gidx = 1 : 2	  % HW/OW
        dc = squeeze(grpmean_DC(gidx,:,steps(mod,:)));
        dc(isinf(dc)|isnan(dc)) = 0;
        obj = plot_hemispheres2(squeeze(dc(1:210,:)), {surf_lh,surf_rh}, 'parcellation', parc);
        colormap(obj.figure, cool_warm)
        saveas(gcf, [outpath, 'figures/dc_', group_name{gidx}, num2str(mod), '.png'])
        close(gcf)
    end
end

% subcortex
for gidx = 1 : 2
    for step = 1 : Nstep
        dc_subcor = mean_subcortical(squeeze(grpmean_DC(gidx, :, step))');
        plot_subcortical(dc_subcor, 'ventricles', 'False', 'cmap', 'RdBu_r', 'color_range', dc_range(step,gidx*2-1:gidx*2));
        enigma_colormap(cool_warm)
        saveas(gcf, [outpath, 'figures/dc_subcor_', group_name{gidx}, num2str(step), '.png'])
        close(gcf)
    end
end


%% 3) Hub regions
load([outpath, 'groupdiff_hub.mat'])

% cortex
for gidx = 1 : 2
    for mode = 1 : 2
        if mode == 1
            hub_4steps = hub{gidx}(1:210, 1:4);
        else
            hub_4steps = hub{gidx}(1:210, 5:7);
        end
        obj = plot_hemispheres2(hub_4steps, {surf_lh,surf_rh}, 'parcellation', parc);
        colormap(obj.figure, cool_warm(129:230, :))
        saveas(gcf, [outpath, 'figures/hub_', group_name{gidx}, num2str(mode), '.png'])
        close(gcf)
    end
end

% subcortex
for gidx = 1 : 2
    for step = 1 : Nstep
        hub_1step = hub{gidx}(:, step);
        plot_subcortical(hub_1step(211:end), 'ventricles', 'False', 'cmap', 'RdBu_r', 'color_range', [-1.3 1.3])
        enigma_colormap(cool_warm)
        saveas(gcf, [outpath, 'figures/hub_subcor_', group_name{gidx}, num2str(step), '.png'])
    end
end


%% 4) Group difference: roi-level
load([outpath, 'groupdiff_ROI_ttest.mat']);
significant_T = H .* T;

% cortex
steps = [[1, 2, 3, 4]; [5, 6, 7, 7]];
clim = [[-5 -5 -5 -5]; [5 5 5 5]];
for mod = 1 : 2
    obj = plot_hemispheres2(significant_T(1:210,steps(mod,:)), {surf_lh,surf_rh}, 'parcellation', parc, 'clim', clim);
    colormap(obj.figure, cool_warm)
    saveas(gcf, [outpath, 'figures/ROIttest', num2str(mod), '.png'])
    close(gcf)
end

% subcortex
for gidx = 1 : 2
    for step = 1 : Nstep
        plot_subcortical(significant_T(211:end, step), 'ventricles', 'False', 'cmap', 'RdBu_r', 'color_range', [-5 5])
        enigma_colormap(cool_warm)
        saveas(gcf, [outpath, 'revision/figures/highlow/ROIttest_subcor_', group_name{gidx}, num2str(step), '.png'])
        close(gcf)
    end
end


%% 5) Group difference: network-level
load([outpath, 'groupdiff_NET_ttest.mat']);
red = cool_warm(129+128/2, :);
blue = cool_warm(128/2, :);

for nidx = 1 : num_network
    clf, figure(1)
    hold on
    for step = 1 : Nstep
        h = bar(step, T(step, nidx));
        if T(step, nidx) > 0
            set(h, 'FaceColor', red);
        else
            set(h, 'FaceColor', blue);
        end
    end
    set(gca,'FontSize',15), ylim([-4, 4]), yticks([-4:1:4]), xticks([1:1:7])
    xlabel('Steps'), ylabel('t-statitic')
    hold off
    saveas(gcf, [outpath, 'figures/NETttest_', num2str(nidx),'.png']);
end
