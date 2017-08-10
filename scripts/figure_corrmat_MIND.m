function h = figure_corrmat_COW(matrix_orig,atlas_params,varargin)
% matrix_orig: a NxN (2D) matrix
% atlas_params: a structure containing various pieces of information about
% your ROIs, including
%   num_rois - the number of ROIs
%   networks - a cell array of network label strings
%   colors - a nets X 3 array of colors for the networks
%   sorti - the sorting order for the ROIs into networks
%   tansitions - the transition boundaries marking one network off from the
%   next
%   centers - the center position for each network (for labeling the name)
% varargin: if provided, upper and lower limits on the color scale
%
% Originally from  T. Laumann
% Edited by C. Gratton to make more flexible
% Edited by C. Gratton for Nebraska, 7.7.17

networks = atlas_params.networks;
colors_new = atlas_params.colors;

matrix = matrix_orig(atlas_params.sorti,atlas_params.sorti);

h = figure('Color',[0.8275 0.8275 0.8275],'Position',[56 143 1295 807]); %[56 143 1095 807]

if nargin>2
    if ~isempty(varargin{1}) & ~isempty(varargin{2})
        climlow = varargin{1};
        climhigh = varargin{2};
        imagesc(matrix,[climlow climhigh]);
    else
        imagesc(matrix);
    end
else
    imagesc(matrix);
end

vline_new(atlas_params.transitions,'k',3);
hline_new(atlas_params.transitions,'k',3);
tickpos = atlas_params.centers;
ax = axis;

set(gca,'XTick',tickpos,'Xlim',[ax(1) ax(2)]);

set(gca,'XTicklabel','');
set(gca,'YTicklabel','');    

tx= text(tickpos,ones(1,length(tickpos))*(atlas_params.num_rois+1),networks);
set(tx,'HorizontalAlignment','right','VerticalAlignment','top','Rotation',45);

for i = 1:length(tx)
     set(tx(i),'Color',[colors_new(i,1) colors_new(i,2) colors_new(i,3)],'FontName','Helvetica','FontSize',10,'FontWeight','bold');   
end
set(gca,'FontWeight','bold','FontSize',10);


ty= text(-1*ones(1,length(tickpos)),tickpos-5,networks);
set(ty,'HorizontalAlignment','right','VerticalAlignment','top')

for i = 1:length(ty)
    set(ty(i),'Color',[colors_new(i,1) colors_new(i,2) colors_new(i,3)],'FontName','Helvetica','FontSize',10,'FontWeight','bold');   
end
colorbar;
set(gca,'FontWeight','bold','FontSize',10);

if nargin>4
    if ~isempty(varargin{4})
        title(titletext,'FontWeight','bold','FontSize',10);
    end
end

axis square;