function [PC,WD,deg,mat] = module_metrics_Dartmouth(mat,assignments,thresh,varargin)
% function [PC,WD] = module_metrics(mat,communities)
%
% Function that calculate the PC and WD values for each region,
% given:
%   mat = a correlation matrix (2D) - must be transformed into a binary
%   matrix for these calculations
%   commmunities = a cell array, w/ one cell per module, w node indices in
%   module
%   thresh = threshold for matrix (assumes kden thresholding)
%   varargin(1) = distance exclusion in mm
%   varargin(2) = ROI file or xyz array
% calculates:
%   PC = participation coefficient
%   WD = within module degree
%   as in Guimera & Ameral, 2005, Nature
%   WD_notz = within module degree, not normalized
%
% CG - 9/8/2014

% First check that matrix is input in 2D:
if length(size(mat)) ~= 2
    error('Input matrix must be 2D');
end

% Prep matrix
% 1. remove close connections to diminish local & motion effects
if length(varargin) > 0
    mat(varargin{2}<varargin{1}) = 0; 
end
% 2. Threshold the matrix (density based)
mat = threshold_the_matrix(mat,thresh);
mat = mat>0; %Binarize the matrix [or switch to weighted computations below...]

% determine communities
communities = make_mod(assignments);

% Some constants/initializations
all_nodes = 1:size(mat,1);
WD = ones(size(mat,1),1)*nan;
PC = ones(size(mat,1),1)*nan;

% quick degree calculation
deg = sum(mat,2); % sum across rows to calculate the degree

% Loop through modules
for c = 1:length(communities)
    
    comm_nodes = communities{c};
    comm_mat = mat(comm_nodes,comm_nodes);

    % WD calculation
    ki  = sum(comm_mat,1); %sum across rows for single # per node (density w/in module connections)
    mean_ks = mean(ki); %average across all nodes in module
    sigma_ks = std(ki); %std across all nodes in module
    zi = (ki - mean_ks)/sigma_ks; %take z-score
    WD(comm_nodes) = zi;
            
    if sum(isnan(zi)) > 0
        if sigma_ks == 0 && ~isnan(mean_ks) %if they're all the same, make them all WD = 0
            WD(comm_nodes) = 0;
        else
            error('nan wd');
        end
    end

end


% PC calculation
k = sum(mat,1); % get total degree for each node
for n = 1:length(all_nodes)
	for c = 1:length(communities)
        comm_nodes = setdiff(communities{c},all_nodes(n)); % all nodes in the community that are not node n
        node_to_comm = mat(all_nodes(n),comm_nodes); %here don't use utri version
        kis(c) = sum(node_to_comm);
    end
    PC(n) = 1 - sum((kis/k(n)).^2);
    
    % for disconnected nodes, set PC = 0
    if isnan(PC(n))
        if k(n) == 0
            PC(n) = 0;
        else
            error('Nan! But not disconnected');
        end
    end
end

end

function mods = make_mod(assign)

nets = unique(assign);

for i = 1:length(nets)
    mods{nets(i)} = find(assign == nets(i));
end

end
