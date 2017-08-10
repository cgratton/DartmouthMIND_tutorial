function mod_metrics_wrapper(atlas)
%%% Quick script for runnign module metrics over a range of thresholds
% Expects to be in corrfile directory

homedir = '/data/cn5/caterina/TaskConn_Methods/all_data/';
if strcmp(atlas,'ParcelCenters')
    restdir = [homedir 'rest_fc_scrubbed_concat_Power/'];
else
    restdir = [homedir 'rest_fc_scrubbed_concat_' atlas '/'];
end
outdir = [restdir 'mod_metrics_' atlas '/'];
if ~exist(outdir,'dir')
    mkdir(outdir);
end

% constants
atlas_params = atlas_parameters(atlas,homedir);
[rois x y z] = textread([homedir atlas_params.atlas_file '.txt'],'%s%f%f%f');
thresh_vals = [0.02:.01:.10];
xd = 20;

% Load correlation matrix and infomap information
switch atlas
    case 'Power'
        mat_fname = 'corrfile_BigBrain264TimOrder_roilist_RMAT_fast_filtertmask_orig.mat';
        infomap_dir = 'corrfile_BigBrain264TimOrder_roilist_Tk002to01in001_S1to29_xd20_INFMAP/'; %note: edit if ch params!
    case 'Parcels'
        mat_fname = 'corrfile_Parcels_711-2b_roilist_RMAT_fast_filtertmask_orig.mat';
        infomap_dir = 'corrfile_Parcels_711-2b_roilist_Tk002to01in001_S1to29_xd20_INFMAP/'; %note: edit if ch params!
    case 'ParcelCenters'
        mat_fname = 'corrfile_modified_ParcelCenters_roilist_RMAT_slow_filtertmask_orig.mat';
        infomap_dir = 'corrfile_modified_ParcelCenters_filtertmask_orig_Tk002to01in001_S1to29_xd20_INFMAP/'; %note: edit if ch params!
end
load([restdir mat_fname]); %loads this into corrmat variable
ave_mat = tanh(mean(atanh(corrmat),3)); % average across subjects
ave_mat_half1 = tanh(mean(atanh(corrmat(:,:,1:2:end)),3));
ave_mat_half2 = tanh(mean(atanh(corrmat(:,:,2:2:end)),3));
infomap_comm = dlmread([restdir infomap_dir 'rawassn.txt']);


% also load data from 120 (Power 2013, Tim data) and about ROIs
switch atlas
    case 'Power'
        power_data_fname = '/data/cn5/caterina/TaskConn_Methods/useful_files/noderole_Power2013.xlsx';
        power_data_full = xlsread(power_data_fname);
        pc_col = 8; %column with summed PC information
        rest120_data_pc = power_data_full(:,8);        
        power_sort = get_sorted_ROI_order(power_data_full(:,2:4),[x y z]); %match col's w 711 space
        rest120_data_pc = rest120_data_pc(power_sort);
        save([outdir 'Rest120_data_pc_sum.mat'],'rest120_data_pc');                
    case {'Parcels'}%,'ParcelCenters'}
        parcels120_data_fname = '/data/cn5/caterina/Atlases/Evan_parcellation/PC/participation_coef_Parcels.mat';
        parcels120_data_full = load(parcels120_data_fname);
        orig_thresholds = [.01:.005:.2];
        %parcels120_data_pc = resample_data(parcels120_data_full.out_data,orig_thresholds,thresh_vals);
        parcels120_data_pc = resample_data2(parcels120_data_full.out_data,orig_thresholds,thresh_vals);
        parcels_sort = atlas_params.sorti;
        parcels120_data_pc_sorted = sum(parcels120_data_pc(atlas_params.sorti,:),2);
        rest120_data_pc = sum(parcels120_data_pc,2);
        save([outdir 'Rest120_data_pc_sum.mat'],'rest120_data_pc');
    case 'ParcelCenters'
        parcels120_data_fname = '/data/cn5/caterina/TaskConn_Methods/all_data/120_proc/mod_metrics_ParcelCenters/mod_metrics.mat';
        parcels120_data_full = load(parcels120_data_fname);
        parcels_sort = atlas_params.sorti;
        parcels120_data_pc_sorted = parcels120_data_full.sum_pc(atlas_params.sorti,:);
        rest120_data_pc = parcels120_data_full.sum_pc;
        save([outdir 'Rest120_data_pc_sum.mat'],'rest120_data_pc');        
end



%%%
% Do calculations
pc = ones(size(ave_mat,1),length(thresh_vals))*nan;
wd = ones(size(ave_mat,1),length(thresh_vals))*nan;
wd_notz = ones(size(ave_mat,1),length(thresh_vals))*nan;

pc_sub = ones(size(corrmat,3),size(ave_mat,1),length(thresh_vals))*nan;
wd_sub = ones(size(corrmat,3),size(ave_mat,1),length(thresh_vals))*nan;
wd_notz_sub = ones(size(corrmat,3),size(ave_mat,1),length(thresh_vals))*nan;

pc_half1 = ones(size(ave_mat,1),length(thresh_vals))*nan; pc_half2 = ones(size(ave_mat,1),length(thresh_vals))*nan;
wd_half1 = ones(size(ave_mat,1),length(thresh_vals))*nan; wd_half2 = ones(size(ave_mat,1),length(thresh_vals))*nan;
wd_notz_half1 = ones(size(ave_mat,1),length(thresh_vals))*nan; wd_notz_half2 = ones(size(ave_mat,1),length(thresh_vals))*nan;

for t = 1:length(thresh_vals)
    comm_mods = make_mod(infomap_comm(:,t));
    [pc(:,t), wd(:,t), wd_notz(:,t)] = module_metrics(ave_mat,comm_mods,thresh_vals(t),xd,atlas_params.roi_file);
    for s = 1:size(corrmat,3)
        [pc_sub(s,:,t) wd_sub(s,:,t) wd_notz_sub(s,:,t)] = module_metrics(corrmat(:,:,s),comm_mods,thresh_vals(t),xd,atlas_params.roi_file);
    end            
    [pc_half1(:,t), wd_half1(:,t), wd_notz_half1(:,t)] = module_metrics(ave_mat_half1,comm_mods,thresh_vals(t),xd,atlas_params.roi_file);
    [pc_half2(:,t), wd_half2(:,t), wd_notz_half2(:,t)] = module_metrics(ave_mat_half2,comm_mods,thresh_vals(t),xd,atlas_params.roi_file);
end
sum_pc = sum(pc,2); sum_wd = sum(wd,2); sum_wd_notz = sum(wd_notz,2);
sum_pc_sub = sum(pc_sub,3); sum_wd_sub = sum(wd_sub,3); sum_wd_notz_sub = sum(wd_notz_sub,3);
sum_pc_half1 = sum(pc_half1,2); sum_wd_half1 = sum(wd_half1,2); sum_wd_notz_half1 = sum(wd_notz_half1,2);
sum_pc_half2 = sum(pc_half2,2); sum_wd_half2 = sum(wd_half2,2); sum_wd_notz_half2 = sum(wd_notz_half2,2);
mean_pc = nanmean(pc,2); mean_wd = nanmean(wd,2); mean_wd_notz = nanmean(wd_notz,2);
fout = [outdir 'mod_metrics.mat'];
save(fout,'sum_pc','sum_wd','sum_wd_notz','sum_pc_sub','sum_wd_sub','sum_wd_notz_sub');

%%%

% Make some quick plots of output
figure('Position',[56 143 1095 807]);
subplot(1,3,1);
plot_per_node(pc,'pc',atlas_params,1);

subplot(1,3,2)
plot_per_node(wd,'wd',atlas_params,0);

subplot(1,3,3)
plot_per_node(wd_notz,'wd (not z)',atlas_params,0);
save_fig(gcf,[outdir 'module_metrics_bynode.png']);

figure;
for m = 1:length(atlas_params.mods)
    plot(mean_pc(atlas_params.mods{m}),mean_wd(atlas_params.mods{m}),'o',...
        'MarkerFaceColor',[atlas_params.colors(m,1) atlas_params.colors(m,2) atlas_params.colors(m,3)],...
        'MarkerEdgeColor',[atlas_params.colors(m,1) atlas_params.colors(m,2) atlas_params.colors(m,3)]); hold on;
end
xlabel('Mean PC across density')
ylabel('Mean WD across density')
title('PC by WD')
hline_new(0,'k-',1);
hline_new(3,'k--',1);
save_fig(gcf,[outdir 'module_metrics_PCxWD.png']);

% and a plot averaging over modules
figure;
subplot(2,1,1);
plot_per_mod(sum_pc,atlas_params.mods,atlas_params.colors,atlas_params.networks,'pc');

subplot(2,1,2);
plot_per_mod(sum_wd_notz,atlas_params.mods,atlas_params.colors,atlas_params.networks,'wd (not z)');
save_fig(gcf,[outdir 'module_metrics_bymod.png']);

% and a plot comparing these results to Power 2013/120 Parcels
figure;
sum_pc = sum(pc,2);
[r,p] = corr(rest120_data_pc,sum_pc);
for m = 1:length(atlas_params.mods)
    plot(rest120_data_pc(atlas_params.mods{m}),sum_pc(atlas_params.mods{m}),'o',...
        'MarkerFaceColor',[atlas_params.colors(m,1) atlas_params.colors(m,2) atlas_params.colors(m,3)],...
        'MarkerEdgeColor',[atlas_params.colors(m,1) atlas_params.colors(m,2) atlas_params.colors(m,3)]); hold on;
end
%plot(power_data_pc,sum(pc,2),'ko');
plot([0 6],[0 6],'k');
xlabel('Rest 120 computed summed PC');
ylabel('Rest summed PC');
title(sprintf('r=%.02f, p=%.02f',r,p));
axis square;
save_fig(gcf,[outdir 'module_metrics_compto120.png']);
close('all');

% and compare the results from single subjects to the mean
sum_pc_sub_mean = mean(sum_pc_sub,1)'; sum_pc_sub_se = std(sum_pc_sub,1)'/sqrt(size(sum_pc_sub,1));
figure('Position',[200 200 1200 400]);
subplot(1,3,1)
sum_pc = sum(pc,2);
[r,p] = corr(sum_pc,sum_pc_sub_mean);
for m = 1:length(atlas_params.mods)
    errorbar(sum_pc(atlas_params.mods{m}),sum_pc_sub_mean(atlas_params.mods{m}), sum_pc_sub_se(atlas_params.mods{m}),...
        'o', 'Color',[atlas_params.colors(m,1) atlas_params.colors(m,2) atlas_params.colors(m,3)]); hold on;
end
%plot(power_data_pc,sum(pc,2),'ko');
plot([0 7],[0 7],'k');
xlabel('Mean matrix summed PC');
ylabel('Ind Sub summed PC');
title(sprintf('r=%.02f, p=%.02f',r,p));
xlim([0 7]); ylim([0 7]);
axis square;

subplot(1,3,2)
sum_pc_sub_mean1 = mean(sum_pc_sub(1:2:end,:),1)'; 
sum_pc_sub_mean2 = mean(sum_pc_sub(2:2:end,:),1)'; 
[r,p] = corr(sum_pc_sub_mean1,sum_pc_sub_mean2);
for m = 1:length(atlas_params.mods)
    plot(sum_pc_sub_mean1(atlas_params.mods{m}),sum_pc_sub_mean2(atlas_params.mods{m}),...
         'o', 'MarkerFaceColor',[atlas_params.colors(m,1) atlas_params.colors(m,2) atlas_params.colors(m,3)],...
        'MarkerEdgeColor',[atlas_params.colors(m,1) atlas_params.colors(m,2) atlas_params.colors(m,3)]); hold on;
end
plot([0 7],[0 7],'k');
xlabel('Half 1 of subjects');
ylabel('Half 2 of subjects');
xlim([0 7]); ylim([0 7]);
axis square;
title(sprintf('Individual Subjects, Half groups r=%.02f',r));

subplot(1,3,3)
[r,p] = corr(sum_pc_half1,sum_pc_half2);
for m = 1:length(atlas_params.mods)
    plot(sum_pc_half1(atlas_params.mods{m}),sum_pc_half2(atlas_params.mods{m}),...
         'o', 'MarkerFaceColor',[atlas_params.colors(m,1) atlas_params.colors(m,2) atlas_params.colors(m,3)],...
        'MarkerEdgeColor',[atlas_params.colors(m,1) atlas_params.colors(m,2) atlas_params.colors(m,3)]); hold on;
end
plot([0 7],[0 7],'k');
xlabel('Half 1 of subjects');
ylabel('Half 2 of subjects');
xlim([0 7]); ylim([0 7]);
axis square;
title(sprintf('Mean matrix, Half groups r=%.02f',r));

save_fig(gcf,[outdir 'module_metrics_compMeantoSub.png']);
close('all');

end


function plot_per_mod(metric,mods,colors,networks,metric_name)

for m = 1:length(mods)
    metric_mod = nanmean(metric(mods{m}));
    metric_mod_se = nanstd(metric(mods{m}))/sqrt(length(mods{m}));
    b = bar(m,metric_mod); 
    set(b,'facecolor',[colors(m,1) colors(m,2) colors(m,3)],'edgecolor',[colors(m,1) colors(m,2) colors(m,3)]); hold on;
    c = errorbar(m,metric_mod,metric_mod_se,'.');
    set(c(1),'color',[colors(m,1) colors(m,2) colors(m,3)]);
    set(c(1),'LineStyle','none');
end

set(gca,'XTicklabel','','Xtick',1:length(networks));
tx= text(1:length(mods),ones(1,length(mods))*0,networks);
set(tx,'HorizontalAlignment','right','VerticalAlignment','top','Rotation',45);
ylabel(['Avg summed ' metric_name ' per module']);

end

function plot_per_node(metric,metric_name,atlas_params,plot_labels)

networks = atlas_params.networks;
colors = atlas_params.colors;
mods = atlas_params.mods;
sorti = atlas_params.sorti;
transitions = atlas_params.transitions;
centers = atlas_params.centers;

imagesc(metric(sorti,:))
xlabel('density');
ylabel('ROI');
colorbar()
ax = axis;

% transitions = [1];
% for m = 1:length(mods)
%     transitions = [transitions mods{m}(end)];
% end
% centers = transitions(1:end-1) + ((transitions(2:end) - transitions(1:end-1))/2);

set(gca,'YTicklabel','') 
if plot_labels
    set(gca,'Ytick',centers)
    ty= text(-1*ones(1,length(centers)),centers-5,networks);
    set(ty,'HorizontalAlignment','right','VerticalAlignment','top')
    for i = 1:length(ty)
        set(ty(i),'Color',[colors(i,1) colors(i,2) colors(i,3)],'FontName','Helvetica','FontSize',12,'FontWeight','bold');   
    end
    set(gca,'FontSize',12)
end
hline_new(transitions+.5,'w',3)

title(metric_name)

end

function data_resamp = resample_data(data,start_thresholds,end_thresh_vals)

samp = end_thresh_vals(2) - end_thresh_vals(1);
for e = 1:length(end_thresh_vals)
    th_low = end_thresh_vals(e) - samp/2;
    th_high = end_thresh_vals(e) + samp/2;
    bin = (start_thresholds > th_low) .* (start_thresholds < th_high);
    data_resamp(:,e) = mean(data(:,logical(bin)),2);
end
    
end

function data_resamp = resample_data2(data,start_thresholds,end_thresh_vals)
% This function assumes that there is a single start_threshold
% corresponding to each end_thresh_vals

for e = 1:length(end_thresh_vals)
    tind = find(abs(start_thresholds - end_thresh_vals(e)) < 0.000001); %to deal w weird behavior
    %tind = find(start_thresholds == end_thresh_vals(e)); %to deal w weird behavior
    if length(tind) < 1
       error('No corresponding threshold');
    elseif length(tind) > 1
        error('Too many corresponding thresholds');
    end
    data_resamp(:,e) = data(:,tind);
end


end


function mods = make_mod(comm_assign)

ids = unique(comm_assign);
for i = 1:length(ids)
    mods{i} = find(comm_assign == ids(i));
end

end

function consensusmap = assgn_figure(net_assn,atlas_params,colors,thresholds)

mincol = 2; % column to use for main assignments

figure('Position',[0 0 800 1000]);

subplot(1,4,1:3)
imagesc(net_assn);
colormap(colors);
%axis square;

% add in labels
set(gca,'YTicklabel','')
set(gca,'Ytick',atlas_params.centers)
ty= text(-1*ones(1,length(atlas_params.centers)),atlas_params.centers-2,atlas_params.networks);
set(ty,'HorizontalAlignment','left','VerticalAlignment','top')
for i = 1:length(ty)
    set(ty(i),'Color',[atlas_params.colors(i,1) atlas_params.colors(i,2) atlas_params.colors(i,3)],'FontName','Helvetica','FontSize',8,'FontWeight','bold');
end
set(gca,'FontSize',10)
hline_new(atlas_params.transitions+.5,'k',3)

% and some axis labels
xlabel('Density')
set(gca,'XTicklabel',thresholds);


%%% Calculate consensus
consensusmap = net_assn(:,mincol);
unassigned = find(consensusmap<2); %1 is unassigned
for unassignedindex = unassigned'
    thisassignments = net_assn(unassignedindex,mincol:end);
    thisassignments(thisassignments<2) = [];
    if ~isempty(thisassignments) && length(unique(thisassignments))<3 &&  length(thisassignments)>2
        % enfore that we only assign something new that is reasonably consistent, and
        % is present at more than 2 thresholds
        consensusmap(unassignedindex) = thisassignments(1);
    end
end

subplot(1,4,4)
imagesc(consensusmap);
colormap(colors);
set(gca,'YTicklabel','');
set(gca,'XTicklabel','');
hline_new(atlas_params.transitions+.5,'k',3);
title(['Consensus, based on col=' num2str(mincol)]);

end
function net_colors = get_colors(infomap_comm,assign_type)

switch assign_type
    case 'random'
        net_names = unique(infomap_comm);
        net_colors = distinguishable_colors(length(net_names) - 1);
        net_colors = [1 1 1; net_colors]; % make 1 (unlabeled) white
        
    case 'preassigned'
        net_colors = [1 1 1;...
            1 0 0;...
            0 0 1;...
            0 1 1;...
            .5 0 .5;...
            0 0 0;...
            1 0 1;...
            .9 .9 0;...
            0.7500 0.2500 0;...
            0.5862 0.8276 0.3103;...
            0 1 0;...
            1 .5 0;...
            0.5000 0.5000 0.5000;...
            0.5172 0.5172 1.0000;...
            0 0.5000 0.5000;...
            0.75 0 0.25;...
            0.9655 0.0690 0.3793;...
            1.0000 0.7586 0.5172;...
            0 0 0.4828;...
            0.1379 0.1379 0.0345]; %...
            %0 0.3448 0];%0 0 0.4828; 0.6207 0.3103 0.2759
        
        figure;
        for i = 1:size(net_colors,1)
            aa(1,1,1:3) = net_colors(i,:);
            subplot(5,5,i);
            image(aa);
            title(i);
        end
end

end