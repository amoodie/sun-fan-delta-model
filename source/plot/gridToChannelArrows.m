function gridToChannelArrows(grid, varargin)
% gridToChannelArrows  show the current channel network
    % shows the network as a series of small red arrows.
    %
    % inputs:
    %   grid
    %      A grid object, or a string pointing to a grid object in a .mat
    %      file on disk
    %

    % process a string (file on disk) to a grid
    if ischar(grid)
        env = load(grid);
        grid = env.grid;
        figure(); hold on;
        imagesc(grid.z)
    end

    % handle plotting with various widths
    if ~isempty(varargin)
        error('Not Implemented.')
        % this is a array to use to set widths of lines
        % array = varargin{1};
    end
    
    % a helper array to track where has been visited
    %   this is needed to avoid plotting looped networks
    visitedCells = false(grid.size);
    
    % what is the starting point of the network?
    %   todo: flexibility for multiple inputs?
    channelStartIndices = grid.inletCell(1); % start the walk from the inlet
    [iStart, jStart] = ind2sub(grid.size, channelStartIndices);

    % start the walking and plotting algo from the inlet
    [~] = walkChannelToNodeDrawArrows(grid, iStart, jStart, visitedCells);
        
end

function [visitedCells] = walkChannelToNodeDrawArrows(grid, iStart, jStart, visitedCells)
% walkChannelToNodeDrawArrows  walk the channel pathway and plot
    % walks down a channel pathway and plots arrows. When a branch is
    % encountered, the algo is called recursively. In this way, the whole
    % network is plotted.

    % convenience
    gridsize = grid.size;
    i = iStart;
    j = jStart;

    % walk until find branch or outlet (==1)
    takeStep = true;
    while takeStep

        if visitedCells(i,j)
            % if we have been here before, do nothing and kill the loop
            %    this prevents from plotting looped networks
            takeStep = false;

        else
            % we have not been here before

            % get where to flow to
            ijFlowsTo = grid.nghbrs(grid.flowsToGraph(:, i, j), i, j);

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
                    [visitedCells] = walkChannelToNodeDrawArrows(grid, x, y, visitedCells);
                end
                takeStep = false; % no more walking here
            elseif numel(ijFlowsTo) == 0
                
                takeStep = false; % no more walking here
            end
        end
    end
end
