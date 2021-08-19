function gridLabelChannelInds(grid)

    for k=1:numel(grid.x)
        if grid.channelFlag(k)
            [x,y] = ind2sub(grid.size, k);
            text(y,x,num2str(k),'FontSize',8);
        end
    end
end