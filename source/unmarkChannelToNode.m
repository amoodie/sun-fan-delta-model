function [grid] = unmarkChannelToNode(grid, startIndex, prevIndex, abandonAll)
%
%   startIndex =  the index of the head cell of this pathway
%   prevIndex  =  the index of the cell that flowed to startIndex
%   abandonAll =  flag of whether to abandon everything downstream, or stop
%                 when another pathway is found

    % preallocate the objects to be returned
    newStarts = []; % where to go
    newPrevs = []; % where came from

    % now walk the channel
    walk = true;
    currentIndex = startIndex; % where along the path we currently are
    while walk
        %% see what kind of cell currentIndex is, w.r.t. upstream
        % check what kind of cell currentIndex is with respect to the
        % upstream direction. Look at the flowsFrom array. This tells us
        % whether the cell has flow coming in from somewhere else (is it a
        % confluence or conduit?).
        numberFlowsFrom = sum(grid.flowsFromGraph(:, currentIndex));
        if numberFlowsFrom >= 2
            % this is a channel cell with flow coming into it from two
            % places. Thus, it will still have discharge even after it no
            % longer receives any flow from this abandoned pathway. So, we
            % remove the abandoned `flowsFrom` index, and leave the while
            % loop without adding anything else to walk (because anywhere
            % else we flow to will has discharge (i.e., is not abandoned)).
            disconnectIndex = grid.nghbrs(:, currentIndex) == prevIndex;  % which connection to disconnect
            grid.flowsFromGraph(disconnectIndex, currentIndex) = 0; % set the disconnectIndex cell to 0 to disconnect it

            % check the status of abandonAll to determine whether to stop
            % or continue abandoning along the channel pathway
            if abandonAll
                % unmark the channel state of the branching cell
                grid.channelFlag(currentIndex) = false;
            else
                walk = false; % no more walking, stopped by break below
                break % force the while to break without continuing
            end

        elseif numberFlowsFrom <= 1
            % this is a channel cell with flow from *one or zero* places,
            % so we want to disconnect it from `flowsFrom` index, and then
            % proceed to the next cell(s) in `flowsTo`.

            % disconnect the flowsFrom
            grid.flowsFrom{currentIndex} = []; % reset cell to empty list
            grid.flowsFromGraph(:, currentIndex) = 0;

        else
            error('error found in length of flowsFrom indices')
        end

        %% see what kind of cell currentIndex is, w.r.t. downstream
        % check what kind of cell currentIndex is with respect to the
        % downstream direction. Look at the flowsTo array to find out if
        % this is a branch or conduit or outlet
        numberFlowsTo = sum(grid.flowsToGraph(:, currentIndex));
        if numberFlowsTo == 2
            % this cell is a branching cell, so we have to treat
            % specially

            % unmark the channel state of the branching cell
            grid.channelFlag(currentIndex) = false;

            % we will need to walk both pathways downstream the branch,
            % so we need to return a level up (i.e., leave this
            % function) to restart the walking a pathway process.

            % get the heads of the `newStart` locations for future
            % pathway walks, also need to duplicate the current
            % location (i.e., the branch) as the `newPrevs` of where
            % the `newStarts` received flow from.
            startsBool = grid.flowsToGraph(:, currentIndex);  % true where flowsTo has channels
            newStarts = grid.nghbrs(startsBool, currentIndex);  % the cell indices of the next cells
            
            %newStarts = grid.flowsTo{currentIndex};
            newPrevs = [currentIndex; currentIndex];

            % clear info on where this node would go.
            % (we already recorded where it would go as newStarts)
            grid.flowsToGraph(:, currentIndex) = 0;
            grid.flowsToFrac(:, currentIndex) = 0;

            for i=1:2
                % walk each branch, recursively
                [grid] = unmarkChannelToNode(grid, newStarts(i), newPrevs(i), abandonAll);
            end
            walk = false; % no more walking here

        elseif numberFlowsTo == 1
            % this is a conduit channel (one (or zero) inflow, one
            % outflow), so we just disconnect and then continue down
            % the pathway

            % grab the next step
            nextBool = grid.flowsToGraph(:, currentIndex);  % true where flowsTo goes next
            nextIndex = grid.nghbrs(nextBool, currentIndex);  % the cell indices of the next cells

            % now unset the flowsTo, channelFlag, partition
            grid.flowsToGraph(:, currentIndex) = 0;
            grid.channelFlag(currentIndex) = false;
            grid.flowsToFrac(:, currentIndex) = 0;

            % update where the flow came from to where we are, and then
            % choose the next index as the current index (i.e., take
            % the step down the pathway).
            prevIndex = currentIndex;
            currentIndex = nextIndex;

        elseif numberFlowsTo == 0
            % this is a channel outlet, there is nowhere else to flow
            % to. Unmark the channel and then break the loop (so to
            % leave the function).
            grid.channelFlag(currentIndex) = false;
            walk = false; % break the loop

        else
            error('error found in number of flowsTo indices.')

        end

    end
end
