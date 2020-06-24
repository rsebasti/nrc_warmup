function ax=define_position(rw,cl)
%plotting
rows=rw;cols=cl;
% define spacing parameters
left_margin = 0.05; % left spacing on the figure
bottom_margin = 0.08; % bottom spacing on the figure
right_margin = 0.04; % right spacing on the figure
top_margin = 0.08; % top spacing on the figure
ax_spacing_horizontal = 0.02; % spacing between columns of axes
ax_spacing_vertical = 0.01; % spacing between rows of axes
ax_width = (1-left_margin-right_margin-(cols-1)*ax_spacing_vertical)/cols; % width of each axis
ax_height = (1-top_margin-bottom_margin-(rows-1)*ax_spacing_horizontal)/rows; % height of each axis
[r,c] = ndgrid(1:rows,1:cols);
ax=arrayfun(@(r,c)axes('Position',[left_margin+(c-1)*(ax_width+ax_spacing_vertical) bottom_margin+(r-1)*(ax_height+ax_spacing_horizontal) ax_width ax_height]),r(:),c(:),'UniformOutput',true);