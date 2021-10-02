function subcor = mean_subcortical(bna)
% 
% This function take average feature values of 14 subcortical regions from 36 subcortical regions in BNA atlas.
% (amygdala, hippocampus, globus pallidus, nucleus accumbens, putamen, caudate, and thalamus for left/right hemisphere)
% 
% Inputs :  bna    - Feature vector/matrix of whole BNA regions (246 X N)
% Outputs:  subcor - Feature vector/matrix of 14 subcortical subregions (14 X N)
% 

numval = size(bna, 2);
subcor = zeros(14, numval);
% L
subcor(1,:) = bna(223,:);                         % L-accumbens
subcor(2,:) = mean([bna(211,:); bna(213,:)]);     % L-amygdala
subcor(3,:) = mean([bna(219,:); bna(227,:)]);     % L-caudate
subcor(4,:) = mean([bna(215,:); bna(217,:)]);     % L-hippocampus
subcor(5,:) = bna(221,:);                         % L-pallidum
subcor(6,:) = mean([bna(225,:); bna(229,:)]);     % L-putamen
subcor(7,:) = mean([bna(231,:); bna(233,:); bna(235,:); bna(237,:); bna(239,:); bna(241,:); bna(243,:); bna(245,:)]);	% L-thalamus
% R
subcor(8,:) = bna(224,:);                         % R-accumbens
subcor(9,:) = mean([bna(212,:); bna(214,:)]);     % R-amygdala
subcor(10,:) = mean([bna(220,:); bna(228,:)]);	  % R-caudate
subcor(11,:) = mean([bna(216,:); bna(218,:)]);	  % R-hippocampus
subcor(12,:) = bna(222,:);                        % R-pallidum
subcor(13,:) = mean([bna(226,:); bna(230,:)]);    % R-putamen
subcor(14,:) = mean([bna(232,:); bna(234,:); bna(236,:); bna(238,:); bna(240,:); bna(242,:); bna(244,:); bna(246,:)]);	% R-thalamus
end
