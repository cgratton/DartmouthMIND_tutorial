function hub_colors = hub_colormap(hubvals)

hub_colors = zeros([size(hubvals,1),3]);
map_vals = colormap(jet(101)); close('all');
map_range = [min(hubvals(:)) max(hubvals(:))];
for r = 1:size(hub_colors,1)
    hub_colors(r,:) = map_vals(round((hubvals(r)-map_range(1))/(map_range(2) - map_range(1))*100)+1,:); %add 1 to deal with 0 indexing
end

end