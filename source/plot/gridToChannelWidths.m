function gridToChannelWidths(grid, varargin)

    if ischar(grid)
        env = load(grid);
        grid = env.grid;
    end

    if length(varargin) > 0
        % this is a array to use to set widths of lines
        array = varargin{1};
    end
    
    channelCellList = {};
    channelStartIndices = [4950]; % start the walk from the inlet
    
    % while there are any indices remaining
    while numel(channelStartIndices) > 0

%         channelStartIndex = channelStartIndices(1);
        [iStart, jStart] = ind2sub(grid.size, channelStartIndices(1));
         % channelStartIndex(1);
        %jStart = channelStartIndex(2);
        channelStartIndices = channelStartIndices(2:end); % clip the front off
        
        [channelXY, newStarts] = walkChannelToNode(grid.flowsTo, iStart, jStart);
        channelCellList = [channelCellList; channelXY];
        
        channelStartIndices = [channelStartIndices; newStarts];
        
    end
    
    figure(); hold on;
    imagesc(grid.z)
    for i=1:length(channelCellList)
        channel = channelCellList{i};
        plot(channel(:,2), channel(:,1), '-')
%         xlim([0.5, 10])
%         ylim([0.5, 10])
    end
        
end

function [channelXY, ijFlowsTo] = walkChannelToNode(flowsTo, iStart, jStart)

    gridsize = size(flowsTo);

    i = iStart;
    j = jStart;
    channelXY = [i, j]; % initialize with the start point
    
    takeStep = true;
    
    % walk until find branch or outlet (==1)
    while takeStep
        % get where to flow to
        ijFlowsTo = flowsTo{i, j};
        
        % check if should take another step
        if numel(ijFlowsTo) ~= 1
            takeStep = false;
            break
        end
        
        % take that step
        [i, j] = ind2sub(gridsize, ijFlowsTo);
        
        % add the next step to the array
        channelXY = [channelXY; [i, j]];
        

    end

end