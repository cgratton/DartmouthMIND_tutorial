function [h,n] = make_spring_fig_MIND(corrmat,mat_thresh,colors,varargin)
% function make_spring_fig(corrmat,mat_thresh,atlas_params)
%
% input:
% corrmat: a 2D ROIxROI matrix
% mat_thresh: the density at which to threshold the matrix (e.g., 0.1 = top
% 10% of edges kept)
% atlas_params: a mat file with information about the ROIs and networks for
% coloring the spring embedded plot
%
% Based on scripts originally from T. Laumann
% Edited by C. Gratton


% Spring embedding constants
% Can play with these if you like
kc = 1; % spring constant
L0 = 25; % circle radius


%transform matrix back into regular correlation values for spring embedding
%func (to make it easier to transform into distances...)
if max(corrmat(:)) > 1
    corrmat = tanh(corrmat);
end

%do thresholding by density, not set value!
%this will also threshold diagonal
adj_mat = threshold_the_matrix(corrmat,mat_thresh);

% do the spring embedding
[h,n] = spring_embedding_func_MIND(double(adj_mat),colors,kc,L0);
axis('off');

%add text to the top of the plot with information about each network
if length(varargin) > 0
    for m = 1:length(varargin{1})
        text(-1.2*L0,L0-m*L0/10,['\color[rgb]{' num2str(varargin{2}(m,:)) '}' varargin{1}{m}]);
    end
end

end