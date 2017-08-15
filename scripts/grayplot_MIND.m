function grayplot_MIND(QCinit,whattomask,stage,subject)
% function grayplot_MIND(QC,whattomask,stage,subject)
%
% This function will make grayplots to look at timeseries after different
% stages of procesing. Additionally allows for bad timepoints to be masked
% out from visualization.
% This will allow one to inspect data for artifacts.
%
% QCfile: QC variable (usually saved under 'QC.mat')
%   this file contains timeseries from graymatter, along with
%   information about movement, processing, etc.
% whattomask: 'tmask' or 'none'
% stage: stage of processing, 1-7
%   1 = original pre-processed data
%   2 = demeaned, detrended
%   3 = residuals from nuisance regression (GM, WM, CSF, motion, + derivatives)
%   4 = interpolation of high-motion censored frames
%   5 = bandpass temporal filter (0.009 - 0.08 Hz)
%   6 = demean and detrend again
%
% subject: subject number, 1-10 [for MSC actually sessions from same
% subject]
%
% Note that this function will make figures and save them to the current
% directory.
% Also, note that current colorscale limits assume mode 1000 normalization
% for timeseries, and show 2% signal change
%
% Each of these conventions could be edited as desired for personal use.
%
% Made by CGratton, 5/21/14
% Edited for MIND, 8/1/17
% Based on FCPROCESS code, v4 (JD Power)

QC = QCinit(subject);
% this needs to have: 
% - gray ts, white ts, and CSF ts (these are subsampled here)
% - original mvm params

% output information:
currDir = pwd;
outDir = [currDir '/output/'];
if ~exist(outDir)
    mkdir(outDir);
end

% constants: 
processes = {'orig','demeantrend1','regress','interp','bpfilt','demeantrend2'};
numpts=size(QC.GMtcs,2); %number of timepoints
rightsignallim = [-20:20:20]; % GS and main plot signal limits - 2% assuming mode 1000 normalization
leftsignallim = [0:10:20];  % DVars limits
rylimz=[min(rightsignallim) max(rightsignallim)]; 
lylimz=[min(leftsignallim) max(leftsignallim)];
FDmult = 10; %multiplier to FD to get in range of DVars values
FDthresh = 0.2; % FD threshold to mark frame for scrubbing

% compute data quality metrics --- CG: COMPUTE YOURSELF
[mvm ddt_mvm FD] = compute_FD(QC.MVM);
DVars = compute_DVARS(QC.GMtcs(:,:,stage));
GS = compute_GS(QC.GMtcs(:,:,stage));

%%% create the figure
figure('Position',[1 1 800 1000]);

% plot individual mvm params
subplot(9,1,1:2);
pointindex=1:numpts;
plot(pointindex,mvm);
xlim([0 numpts]); ylim([-1,1]);
ylabel('mvm-XYZPYR');

% Next, plot FD, DVars, and GS on the same plot
subplot(9,1,3:4);
[h a(1) a(2)]=plotyy(pointindex,DVars,pointindex,GS);
set(h(1),'xlim',[0 numpts],'ylim',lylimz,'ycolor',[0 0 0],'ytick',leftsignallim);
ylabel('R:FD*10 B:DV G:GS');
set(a(1),'color',[0 0 1]);
set(h(2),'xlim',[0 numpts],'ycolor',[0 0 0],'xlim',[0 numel(DVars)],'ylim',rylimz,'ytick',rightsignallim);
set(a(2),'color',[0 .5 0]);
axis(h(1)); hold on;
plot([1:numpts],FD*FDmult,'r');
hline_new(FDthresh.*FDmult,'k',1);
%h4=refline([FDthresh.*FDmult, 0],'k');
hold off;

% next plot gray matter signal, masking if requested
subplot(9,1,5:8);
new_GMtcs = QC.GMtcs(:,:,stage);
if strcmp(whattomask,'tmask')
    new_GMtcs(:,QC.tmask==0) = min(QC.GMtcs(:)); %set all masked values to large negative value so black
end
imagesc(new_GMtcs(:,:),rylimz); colormap(gray); ylabel('GRAY'); %CG CHANGES END here

% finally, plot WM and CSF ts
subplot(9,1,9);
imagesc([QC.WMtcs(:,:,stage);QC.CSFtcs(:,:,stage)],rylimz); ylabel('WM CSF');


% save out
saveas(gcf,[outDir 'Sess' num2str(subject) '_' num2str(stage) '_' processes{stage} '_mask' whattomask '.tiff'],'tiff');


end

function [mvm ddt_mvm FD] = compute_FD(mvm_orig)

% convert rotational mvm params from deg to mm, assuming 50mm head size [has already been
% done in this case]
% radius = 50;
% mvm(:,1:3) = mvm_orig(:,1:3); % translation params stay the same
% mvm(:,4:6) =mvm_orig(:,4:6).*(2*radius*pi/360); % rotation params converted (calculating arc segment length)
mvm = mvm_orig;

% take original movement parameters, demean and detrend
ddt_mvm = diff(mvm);
ddt_mvm = [zeros(1,6); ddt_mvm];

% compute FD
FD=(sum(abs(ddt_mvm),2));

end


function DVARS = compute_DVARS(GMtcs)

DVARS = rms(diff(GMtcs,1,2));
DVARS = [0 DVARS];

end

function GS = compute_GS(GMtcs)

GS = nanmean(GMtcs,1);

end