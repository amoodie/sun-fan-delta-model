function gridLabelChannelInds(grid)
%gridLabelChannelInds Add cell index labels to a plot
% Adds the cell indices of cells marked by `grid.channelFlag` to the
% current axis.
    
    for k=1:numel(grid.x)
        if grid.channelFlag(k)
            [x,y] = ind2sub(grid.size, k);
            text(y,x,num2str(k),'FontSize',8);
        end
    end
end
