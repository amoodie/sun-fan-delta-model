function grid=unmarkAbandonedChannels(grid,Qw_threshold)
%unmarkAbandonedChannels Unmark abandoned channels in tracking arrays
%   This function loops through the `grid`, and along channel pathways, to
%   unmark (i.e., disconnect) channel pathways that have been abandoned.
%   The returned `grid` is updated so that locations that are determined as
%   abandoned (based on Qw_threshold) no longer contain `flowsTo` or
%   `flowsFrom` information and `channelFlag == false`.
%
%   glossary:
%    - branch: a cell that flows to two cells
%    - channel: a cell that is marked as `channelFlag==true`
%    - pathway: a series of cells connected by `flowsTo`
%    - conduit: a cell that has one inflow (one index in `flowsFrom`) and
%                 one outflow (one index in `flowsTo`.

    % make a list of all the branches in the domain
    branchCellIndices = []; % a list of where a channel cell flows to two locations (i.e., branch cells)
    for k=1:numel(grid.channelFlag)
        if numel(grid.flowsTo{k}) == 2
            branchCellIndices = [branchCellIndices, k]; %#ok<*AGROW>
        end
    end

    % loop through the branches and check cells downstream of the branch
    % point (i.e., the `flowsTo{branchIndex}`) against the discharge
    % threshold.
    for i=1:numel(branchCellIndices)

        branchIndex = branchCellIndices(i); % the ith branch index

        % try both pathways against the threshold discharge
        for j=1:numel(grid.flowsToFrac_Qw_distributed{branchIndex})
            testDischarge = grid.flowsToFrac_Qw_distributed{branchIndex}(j); % the discharge fraction in the jth branch
            if testDischarge < Qw_threshold
                % this channel pathway needs to be disconnected and unset
                % along the channel course.
                %
                % while we walk this path, we may encounter other branches,
                % as well as channel cells that receive water from multple
                % sources. To make sure we unmark all abandoned channels in
                % the domain, we walk through the domain and accumulate a
                % list of indices that head pathways that need to be
                % abandoned, and loop through that list until it is empty.
                % Importantly, we also need to track where we entered the
                % current cell *from*, so that we can unset the correct
                % flowsFrom cell when a cell has multiple flowsFrom (only
                % one of which is the abandoned pathway we walked down).
                startIndexList = grid.flowsTo{branchIndex}(j); % cell index list for path starting-point (index of abandoned channel head, not branch)
                prevIndexList = branchIndex; % cell index list for previous cell of path starting-point

                % now, we start to actually do the unsetting and walking.

                % unset flowsTo *at the branch* and update partitioning
                grid.flowsTo{branchIndex} = grid.flowsTo{branchIndex}(1:2 ~= j);
                grid.flowsToFrac_Qw_distributed{branchIndex} = 1;

                % recursively walk the channel path and unset anything
                % that is abandoned (i.e., cells that have no *other*
                % `flowsFrom` than the abandoned path.
                % note: this has to be recursive, because we may encounter branches
                % along the abandoned pathway and we need to clear each of
                % these pathways.
                while numel(startIndexList) > 0

                    % get the new starting cell index and previous cell
                    % index from the list
                    startIndex = startIndexList(1);
                    prevIndex = prevIndexList(1);

                    % update the list by stripping away the first index
                    startIndexList = startIndexList(2:end);
                    prevIndexList = prevIndexList(2:end);

                    % now, recursively walk the pathway from `startIndex`
                    % if a branch is encountered along the walk, then two
                    % downstream cells (newStarts) will be added to the
                    % `startIndexList`, and the branching cell is
                    % duplicated into `newPrevs`.
                    [grid, newStarts, newPrevs] = unmarkChannelToNode(grid, startIndex, prevIndex);
                    startIndexList = [startIndexList; newStarts];
                    prevIndexList = [prevIndexList; newPrevs];

                end % end while loop, all pathways have been walked

                % we removed one pathway (the jth pathway) of this branch
                % (the ith branch), so we don't want to check the other
                % (i.e., skip the other j)
                break % break the `for j=` loop and go to next `i` branch

            end
        end
    end
end

function [grid, newStarts, newPrevs] = unmarkChannelToNode(grid, startIndex, prevIndex)
%
%   startIndex =  the index of the head cell of this pathway
%   prevIndex =  the index of the cell that flowed to startIndex

    % preallocate the objects to be returned
    newStarts = []; % where to go
    newPrevs = []; % where came from

    % now walk the channel
    walk = true;
    currentIndex = startIndex; % where along the path we currently are
    while walk

        if numel(grid.flowsFrom{currentIndex}) >= 2
            % this is a channel cell with flow coming into it from two
            % places. Thus, it will still have discharge even after it no
            % longer receives any flow from this abandoned pathway. So, we
            % remove the abandoned `flowsFrom` index, and leave the while
            % loop without adding anything else to walk (because anywhere
            % else we flow to will has discharge (i.e., is not abandoned)).
            keepidx = (grid.flowsFrom{currentIndex} ~= prevIndex);
            grid.flowsFrom{currentIndex} = grid.flowsFrom{currentIndex}(keepidx);
            walk = false; % no more walking

        elseif numel(grid.flowsFrom{currentIndex}) <= 1
            % this is a channel cell with flow from *one or zero* places,
            % so we want to disconnect it from `flowsFrom` index, and then
            % proceed to the next cell(s) in `flowsTo`.

            % disconnect the flowsFrom
            grid.flowsFrom{currentIndex} = []; % reset cell to empty list

            if numel(grid.flowsTo{currentIndex}) == 2
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
                newStarts = grid.flowsTo{currentIndex};
                newPrevs = [currentIndex; currentIndex];

                % clear info on where this node would go.
                % (we already recorded where it would go as newStarts)
                grid.flowsTo{currentIndex} = []; % reset cell to empty list
                grid.flowsToFrac_Qw_distributed{currentIndex} = [];

                walk = false; % no more walking here

            elseif numel(grid.flowsTo{currentIndex}) == 1
                % this is a conduit channel (one (or zero) inflow, one
                % outflow), so we just disconnect and then continue down
                % the pathway

                % grab the next step
                nextIndex = grid.flowsTo{currentIndex};

                % now unset the flowsTo, channelFlag, partition
                grid.flowsTo{currentIndex} = [];
                grid.channelFlag(currentIndex) = false;
                grid.flowsToFrac_Qw_distributed{currentIndex} = [];

                % update where the flow came from to where we are, and then
                % choose the next index as the current index (i.e., take
                % the step down the pathway).
                prevIndex = currentIndex;
                currentIndex = nextIndex;

            elseif numel(grid.flowsTo{currentIndex}) == 0
                % this is a channel outlet, there is nowhere else to flow
                % to. Unmark the channel and then break the loop (so to
                % leave the function).
                grid.channelFlag(currentIndex) = false;
                walk = false; % break the loop

            else
                error('error found in number of flowsTo indices.')

            end
        else
            error('error found in length of flowsFrom indices')
        end
    end
end
