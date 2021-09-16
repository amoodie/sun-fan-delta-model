function gridToChannelArrows(grid, varargin)

    if ischar(grid)
        env = load(grid);
        grid = env.grid;
        figure(); hold on;
        imagesc(grid.z)
    end

    if ~isempty(varargin)
        % this is a array to use to set widths of lines
        array = varargin{1};
    end
    
    visitedCells = false(grid.size);
    
    channelStartIndices = grid.inletCell; % start the walk from the inlet

    
    [iStart, jStart] = ind2sub(grid.size, channelStartIndices(1));
    [~] = walkChannelToNodeDrawArrows(grid.flowsTo, iStart, jStart, visitedCells);
        
end

function [visitedCells] = walkChannelToNodeDrawArrows(flowsTo, iStart, jStart, visitedCells)

    gridsize = size(flowsTo);

    i = iStart;
    j = jStart;

    takeStep = true;

    % walk until find branch or outlet (==1)
    while takeStep

        if visitedCells(i,j)
            % if we have been here before, do nothing and kill the loop
            takeStep = false;

        else
            % we have not been here before

            % get where to flow to
            ijFlowsTo = flowsTo{i,j};

            p1 = [i,j]; % source point

            % for all places this source flows to (0 or 1 or 2)
            for bb=1:numel(ijFlowsTo)
                % draw this line
                [x,y] = ind2sub(gridsize, ijFlowsTo(bb));
                p2 = [x,y]; % destination point
                dp = p2-p1; % x-y distances between (vector lengths)

                % draw and adjust the arrow
                q = quiver(p1(2),p1(1),dp(2),dp(1),'r','LineWidth',1);
                q.Marker = '.';
            end

            % record that we visited this cell
            visitedCells(i, j) = true;

            % now determine the next step to take
            if numel(ijFlowsTo) == 1
                % take that step
                [i, j] = ind2sub(gridsize, ijFlowsTo);
            elseif numel(ijFlowsTo) == 2
                % a branch
                % traverse any pathways below this branch (recursion)
                for bb=1:2
                    % walk each branch, recursively
                    [x,y] = ind2sub(gridsize, ijFlowsTo(bb));
                    [visitedCells] = walkChannelToNodeDrawArrows(flowsTo, x, y, visitedCells);
                end
                takeStep = false; % no more walking here
            elseif numel(ijFlowsTo) == 0
                
                takeStep = false; % no more walking here
            end
        end
    end
end