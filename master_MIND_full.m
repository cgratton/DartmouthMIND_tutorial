%%% MIND: C Gratton RSFC tutorial

%%% This tutorial looks at different ways of visualizing your data
% 1. First, I present a way developed by JD Power for looking at your data
% at different stages of processing to detect artifacts
% 2. Next, we turn to FC measures. We look at correlation matrices and
% their variability.
% 3. Next we look at graph representations of these correlation matrices
% 4. Finally, we examine hub measures for these graphs.

% Dataset is from Midnight Scan Club, 10 highly sampled individual subjects
% Freely available on OpenfMRI: https://openfmri.org/dataset/ds000224/
% See paper with description:
%   Gordon et al. (2017) Precision Functional Mapping of Individual Human Brains. Neuron
%   http://www.cell.com/neuron/fulltext/S0896-6273(17)30613-X

%% Initialization of directory information:

thisDir = [pwd '/'];
outdir = [thisDir 'output/'];
datadir = [thisDir 'data/'];
scriptdir = [thisDir 'scripts/'];
addpath(scriptdir);

%% Grayplots

% Start with movement parameters and time-series data - opens into
% GrayplotInfo structure
%for now only provided for MSC01, since files can be large.
% See MSC openfMRI dataset for data from additional subjects
load([datadir 'MSC01_grayplot_data.mat']);

% Go over how to calculate FD, DVars, and GS on session 1
% FD: average absolute displacement from frame to frame
% DVARS: root mean sq difference in voxelwise signal from frame to frame
% GS: average signal across all gray matter voxels

% make grayplots for single session, different stages of processing.
% Discuss consequences of each stage of processing for signal.
for st = 1:6
    grayplot_MIND(GrayplotInfo,'none',st,1);
end

% make grayplots for different sessions. Discuss variability.
for i = 1:10
    grayplot_MIND(GrayplotInfo,'none',6,i)
end

% Discuss other finer points:
%   1. How to quantitatively check for motion contamination - QC-FC
%   measures:
%       % Hi vs. Low FD group differences
%       % Correlation between FD and FC (QC-FC)
%       % Distance dependence of FD-FC correlations
%   2. Relatedly, how to set FD threshold
%       % by eye, based on grayplots
%       % comparison (hi vs. low motion groups) between motion censoring
%       and random censoring
%   3. Other issues: 
%       % we also require: 3-5 contiguous frames, 50 frame minimum per run, 120 frames minimum per person
%       %respiration artifact in movement parameters

% Relevant references:
% For more insight into what this sort of approach tells you, read:
%   Power, J. D. (2017). A simple but useful way to assess fMRI scan qualities. Neuroimage, 154, 1, 150-158
% And for more information on FC-processing strategies, read:
%   Power, J.D., et al., (2014). Methods to detect, characterize, and remove motion artifact in resting state fMRI. Neuroimage, 84, 320-341
%   Power J.D., et al. (2015). Recent progress and outstanding issues in motion correction in resting-state fMRI. Neuroimage, 105, 536-551.
%   Power, J.D., et al. (2017). Sources and implications of whole-brain fMRI signals in humans. Neuroimage, 146, 609-625
%   Ciric, R., et al., (2017). Benchmarking of participant-level confound regression strategies for the control of motion artifact in studies of functional connectivity. Neuroimage, 154, 174-187


%% Correlation matrices

load([datadir 'Parcel_params.mat']); % load mat file with information about ROIs and networks

% Start with TS ROIs from 333 and tmask, along with network assignments
% [and possibly also from individual parcels]
%   1. Discuss how we got these ROIs - surface mapping procedure in between
%        [quickly make an image of the ROI timeseries, pre and post mask]
%   2. Discuss different types of ROIs, pros and cons
%   3. Discuss group vs. individual ROIs
for s = 1:10
    ROIdata(s) = load(sprintf('%sMSC%02d_parcel_timecourse.mat',datadir,s));
end

% compute correlations among TS ROIs, masking out bad frames
for s = 1:10
    for i = 1:10
        corrmat(s,i,:,:) = atanh(corr(ROIdata(s).parcel_time{i}(logical(ROIdata(s).tmask_all{i}),:)));
    end
end

% plot a few examples of said correlations - order ROIs randomly
figure;
rand_roi = randperm(Parcel_params.num_rois);
groupmat = squeeze(mean(mean(corrmat,2),1));
imagesc(groupmat(rand_roi,rand_roi),[-0.4,1]);
colormap('jet');
title('Group Matrix, Random ROI order');
colorbar;

% plot correlations in network order for the group
%   1. Discuss/do community detection methods (Infomap, modularity
%   optimization)
%   2. Discuss multi-scale nature of networks
figure_corrmat_MIND(groupmat,Parcel_params,-0.4,1);
title('Group Matrix, ROIs ordered by network');
colormap('jet');
saveas(gcf,[outdir 'Corrmat_group.tiff'],'tiff');
close('all');

% Look at relationship between correlation matrices across sessions and
% subjects
%   1. Look by eye at group and variability
for i = 1:10 %different sessions, same subject
    figure_corrmat_MIND(squeeze(corrmat(1,i,:,:)),Parcel_params,-0.4,1);
    title(['Subject 1, session ' num2str(i)]);
    colormap('jet');
    saveas(gcf,sprintf('%sCorrmat_MSC01_sess%02d.tiff',outdir,i),'tiff');
end
close('all');
for s = 1:10 %different subjects, avg over sessions
    figure_corrmat_MIND(squeeze(mean(corrmat(s,:,:,:),2)),Parcel_params,-0.4,1);
    title(['AllSessAvg, Subject ' num2str(s)]);
    colormap('jet');
    saveas(gcf,sprintf('%sCorrmat_MSC%02d_sessmean.tiff',outdir,s),'tiff');
end
close('all');

% Make similarity matrices
maskmat = ones(Parcel_params.num_rois);
maskmat = logical(triu(maskmat,1));
count = 1;
for s = 1:10
    for i = 1:10
        tmp = corrmat(s,i,:,:);
        corrlin(count,:) = tmp(maskmat);
        count = count+1;
    end
end
simmat = corr(corrlin');
figure('Position',[1 1 1000 800]);
imagesc(simmat,[0 1]); colormap('jet');
hline_new([10:10:90]+0.5,'k',2);
vline_new([10:10:90]+0.5,'k',2);
set(gca,'XTick',[5:10:95],'YTick',[5:10:95],...
    'XTickLabel',{'MSC01','MSC02','MSC03','MSC04','MSC05','MSC06','MSC07','MSC08','MSC09','MSC10'},...
    'YTickLabel',{'MSC01','MSC02','MSC03','MSC04','MSC05','MSC06','MSC07','MSC08','MSC09','MSC10'});
axis square;
colorbar;
title('Correlation Matrix Similarity');
saveas(gcf,[outdir 'SimilarityMat_rest.tiff'],'tiff');
close('all');

% Talk about doing these measures with individual ROIs

% for more information, relevant references:
%   Power et al. (2011). Functional Network Organization of the Human Brain. Neuron, 72, 4, 665-678
%   Infomap: http://www.mapequation.org/code.html
%       Rosvall & Bergstrom (2008). Maps of invofmation flow reveal
%       community structure in complex networks. PNAS, 105, 1118


%% Spring-embedding plots

% Start with a single correlation matrix per subject and network assignment
% [For the group and for the individual separately?]
submats = squeeze(mean(corrmat,2));

% and colormat variable needed for function:
% information about prespecified networks
colors = zeros([size(corrmat,1),3]);
for m = 1:length(Parcel_params.mods)
    colors(Parcel_params.mods{m},:) = repmat(Parcel_params.colors(m,:),[length(Parcel_params.mods{m}),1]);
end

% Create spring embedding plots at one threshold for the group
%   1. Discuss different types of thresholding; discuss weighted vs. not
%   2. Discuss issues of network size and density for feasibility
make_spring_fig_MIND(groupmat,0.02,colors,Parcel_params.networks,Parcel_params.colors);

% Now create across a range of thresholds (1% - 10%, 20%, 30%) for the
% group -- IF TIME
%   1. Discuss consequences of different thresholds
thresholds = [0.01:0.01:0.04]; %[0.01:0.01:0.1,0.2,0.3]; %these take a while, esp for higher thresholds. Only do a few.
for t = 1:length(thresholds)
    %tic;
    make_spring_fig_MIND(groupmat,thresholds(t),colors,Parcel_params.networks,Parcel_params.colors);
    saveas(gcf,sprintf('%sSpring_group_t%.02f.tiff',outdir,thresholds(t)),'tiff');
    %toc
end
close all;

% Pick your favorite threshold, and do the same for all subjects -- IF TIME
%   1. Discuss differences
%   2. Discuss consequences of group vs. individual ROIs and network
%   assignments
for s = 1:10
    make_spring_fig_MIND(squeeze(submats(s,:,:)),0.02,colors,Parcel_params.networks,Parcel_params.colors);
    saveas(gcf,sprintf('%sSpring_sub%02d_thresh0.02.tiff',outdir,s),'tiff');
end
close all;

% for more information, relevant references:
%   Bullmore & Sporns (2009). Complex brain networks: graph theoretical analysis of structural and functional systems. Nature Reviews Neuroscience, 10 (3), 186-198
%   Sporns (2010). Networks of the Brain. MIT Press.
%   Spring Embedding Methods (Kamada-Kawai): https://arxiv.org/pdf/1201.3011.pdf
    


%% Hub measures [If time]

thresholds = [0.01:0.01:0.10];

% Start with a set of thresholded correlation matrices and network assignments for each threshold
%   [In the interest of time/ease, I am pre-computing network assignments
%   and providing them here]
infomapcomm = load([datadir 'Allsubavg_333parcels_infomapassn.mat']);
% code for plotting infomap output:
% colors = distinguishable_colors(max(unique(infomapcomm.clrs)));
% colors(1,:) = [1 1 1];
% figure;
% imagesc(infomapcomm.clrs(Parcel_params.sorti,:),[1 max(unique(infomapcomm.clrs))]);
% hline_new(Parcel_params.transitions,'k',2);
% title('Assignments across thresholds');
% colormap(colors);

% Compute hub measures - degree, PC, and WD - at one threshold in one
% subject
%   1. Look at formula, practice computing measures at one thershold in group
%       degree = total number of edges to a node
%       PC = 1 - for all mods (the # of edges of node to mod m/total # of edges of node)
%       WD = z-score normalized measure of degree within an network
%       See formula from Guimera & Amaral (2005). Functional Cartography of
%       Complex Metabolic Networks. Nature, 433, 895-900
%       http://www.nature.com/nature/journal/v433/n7028/full/nature03288.html?foxtrotcallback=true
%       Look at article end for formulas
%   2. Discuss thresholding, distance exclusion
[pc, wd, degree] = module_metrics_Dartmouth(groupmat,infomapcomm.clrs(:,1),0.02,Parcel_params.dist_thresh,Parcel_params.dmat);

% Now do this across thresholds
%   0. Plot image of hub measures across thresholds with network boundaries
%   1. Discuss how they change across thresholds
%   2. Discuss how they change across networks
%   3. Discuss potential issues with each measure
for t = 1:length(thresholds)
    [pc(:,t), wd(:,t), degree(:,t)] = module_metrics_Dartmouth(groupmat,infomapcomm.clrs(:,t),thresholds(t),Parcel_params.dist_thresh,Parcel_params.dmat);
end
figure_hubs(Parcel_params,thresholds,degree,wd,pc);
saveas(gcf,sprintf('%sHubmeasures_group.tiff',outdir),'tiff');

% [If time, do this across other subjects]
for s = 1:10
    infomapcomm = load(sprintf('data/MSC%02d_333parcels_infomapassn.mat',s));
    for t = 1:length(thresholds)
        [pc_sub(:,t,s), wd_sub(:,t,s), degree_sub(:,t,s)] = module_metrics_Dartmouth(squeeze(submats(s,:,:)),infomapcomm.clrs(:,t),thresholds(t),Parcel_params.dist_thresh,Parcel_params.dmat);
    end
    figure_hubs(Parcel_params,thresholds,degree_sub(:,:,s),wd_sub(:,:,s),pc_sub(:,:,s));
    saveas(gcf,sprintf('%sHubmeasures_MSC%02d.tiff',outdir,s),'tiff');
    clear infomapcomm;
end
close('all');

% [If time] Make a spring embedded plot, colored by hub measures rather than networks
%   1. Start with group and favorite threshold. Do other versions if time.
t = 2;
hub_colors = hub_colormap(pc(:,t));
make_spring_fig_MIND(groupmat,thresholds(t),hub_colors);
saveas(gcf,sprintf('%sSpring_Group_PC_t%.02f.tiff',outdir,thresholds(t)),'tiff');

hub_colors = hub_colormap(wd(:,t));
make_spring_fig_MIND(groupmat,thresholds(t),hub_colors);
saveas(gcf,sprintf('%sSpring_Group_WD_t%.02f.tiff',outdir,thresholds(t)),'tiff');

hub_colors = hub_colormap(degree(:,t));
make_spring_fig_MIND(groupmat,thresholds(t),hub_colors);
saveas(gcf,sprintf('%sSpring_Group_degree_t%.02f.tiff',outdir,thresholds(t)),'tiff');


% And for work on hubs and their importance in brain function, in addition to the references above, see:
%   Gratton, C., et al., (2012). Focal brain lesions to critical locations cause widespread disruption of the modular organization of the brain. Journal of Cognitive Neuroscience, 24 (6), 1275-1285
%   Power, J.D. et al. (2013). Evidence for hubs in human functional brain networks. Neuron, 79 (4), 798-813
%   Warren, D.E., et al. (2014). Network measures predict neuropsychological outcome after brain injury. PNAS, 111 (39), 14247-14252
%
% For a general intro on using graph theory in brain networks see: 
%   Bullmore & Sporns (2009). Complex brain networks: graph theoretical
%       analysis of structural and functional systems. Nature Reviews
%       Neuroscience 10(3), 186.
%   Sporns (2010). Networks of the brain
% For extensions of many graph metrics to weighted networks see also:
%   Rubinov & Sporns, 2011. Weight-conserving characterization of
%       complex functional brain networks. Neuroimage, 56(4) 2068-2079.
%   Rubinov & Sporns, 2010. Complex network measures of brain
%       connectivity: uses and interpretations, 52(3) 1059-1069
%
% The following packages contain tools for graph theoretical analyses:
%   Brain Connectivity Toolbox (Sporns, Matlab/Python/C++): https://sites.google.com/site/bctnet/
%   NetworkX (Python): https://networkx.github.io/ 
%       see also brainx extension: https://github.com/nipy/brainx
