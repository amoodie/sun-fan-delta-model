function [grid] = countFlowsToInds(grid)
%countFlowsToInds Count the number of places each cell flows
    % returns NaN if not a channel, otherwise the count of flowsTo
    
    % everything starts as NaN
    grid.flowsToCount(:) = NaN; 

    % sum down the first axis for all connections
    flowsToCount = squeeze(sum(grid.flowsToGraph, 1));
    
    % only set to count where there are channels (NaN elsewhere)
    grid.flowsToCount(grid.channelFlag) = flowsToCount(grid.channelFlag);

end
