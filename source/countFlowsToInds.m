function [grid] = countFlowsToInds(grid)
%countFlowsToInds Count the number of places each cell flows
    % returns NaN if not a channel, otherwise the count of flowsTo
    
    % everything starts as NaN
    grid.flowsToCount(:) = NaN; 
    
    % loop through channel cells, and determine the number of places it
    % flowsTo
    channelInds = find(grid.channelFlag);
    
    % loop
    for kk=1:numel(channelInds)
        % kth channelInd
        k = channelInds(kk);
        
        % count it
        grid.flowsToCount(k) = numel(grid.flowsTo{k});
    end
    
    flowsToCount = sum(grid.flowsToGraph, 1);
    
    if ~all(flowsToCount(grid.channelFlag) == grid.flowsToCount(grid.channelFlag))
        keyboard
    end
    
%     grid.flowsToCount(grid.channelFlag) = flowsToCount(grid.channelFlag);
    
end

