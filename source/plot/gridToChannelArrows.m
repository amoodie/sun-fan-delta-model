function gridToChannelArrows(grid, varargin)

    if ischar(grid)
        env = load(grid);
        grid = env.grid;
        figure(); hold on;
        imagesc(grid.z)
    end

    if length(varargin) > 0
        % this is a array to use to set widths of lines
        array = varargin{1};
    end
    
    channelCellList = {};
    channelStartIndices = [4950]; % start the walk from the inlet

    % while there are any indices remaining
    while numel(channelStartIndices) > 0

        [iStart, jStart] = ind2sub(grid.size, channelStartIndices(1));
        %jStart = channelStartIndex(2);
        channelStartIndices = channelStartIndices(2:end); % clip the front off
        
        [channelXY, newStarts] = walkChannelToNodeDrawArrows(grid.flowsTo, iStart, jStart);
        channelCellList = [channelCellList; channelXY];
        
        channelStartIndices = [channelStartIndices; newStarts];
        
    end
        
end

function [channelXY, ijFlowsTo] = walkChannelToNodeDrawArrows(flowsTo, iStart, jStart)

    gridsize = size(flowsTo);

    i = iStart;
    j = jStart;
    channelXY = [i, j]; % initialize with the start point
    
    takeStep = true;
    
    % walk until find branch or outlet (==1)
    while takeStep
        % get where to flow to
        ijFlowsTo = flowsTo{i, j};
        
        p1 = [i, j];                         % First Point
        for bb=1:numel(ijFlowsTo)
            [x,y] = ind2sub(gridsize, ijFlowsTo(bb));
            p2 = [x,y];   % Second Point
            dp = p2-p1  ;                       % Difference

            q = quiver(p1(2),p1(1),dp(2),dp(1),'r','LineWidth',1);
            q.Marker = '.';
        end

        % check if should take another step
        if numel(ijFlowsTo) ~= 1
            takeStep = false;
        end
        
        % take that step
        [i, j] = ind2sub(gridsize, ijFlowsTo);
%         
%         % add the next step to the array
%         channelXY = [channelXY; [i, j]];
        

    end

end