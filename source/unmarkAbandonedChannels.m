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
    branchCellIndices = find(grid.flowsToCount >= 2)';

    % loop through the branches and check cells downstream of the branch
    % point (i.e., the `flowsTo{branchIndex}`) against the discharge
    % threshold.
    for i=1:numel(branchCellIndices)

        branchIndex = branchCellIndices(i); % the ith branch index

        % try all pathways against the threshold discharge
        numberFlowsTo = sum(grid.flowsToGraph(:, branchIndex));
        channelsBool = grid.flowsToGraph(:, branchIndex);
        flowsToFracs = grid.flowsToFrac(channelsBool, branchIndex);
        flowsToInds = grid.nghbrs(channelsBool, branchIndex);
        %flowsToIdxs = find(flowsToFracs);
        for j=1:numberFlowsTo
            testDischargeFrac = flowsToFracs(j); % the discharge fraction in the jth branch

            if (testDischargeFrac < Qw_threshold)
                % this channel pathway needs to be disconnected and unset
                % along the channel course.
                %
                % while we walk this path, we may encounter other branches,
                % as well as channel cells that receive water from multple
                % sources. To make sure we unmark all abandoned channels in
                % the domain, we recursively traverse each branch we
                % encounter before returning. Importantly, we also need to
                % track where we entered the current cell *from*, so that
                % we can unset the correct flowsFrom cell when a cell has
                % multiple flowsFrom (only one of which is the abandoned
                % pathway we walked down).                
                startIndex = flowsToInds(j); % cell index for path starting-point (index of abandoned channel head, not branch)
                prevIndex = branchIndex; % cell index for previous cell of path starting-point

                % now, we start to actually do the unsetting and walking.

                % unset flowsTo *at the branch* and update partitioning
                which_index = find((startIndex - branchIndex) == grid.iwalk);
                grid.flowsToGraph(which_index, branchIndex) = 0;
                grid.flowsToFrac(which_index, branchIndex) = 0;
                
                % renormalize the other outputs to equal 1
                renorm = grid.flowsToFrac(:, branchIndex) / sum(grid.flowsToFrac(:, branchIndex));
                grid.flowsToFrac(:, branchIndex) = renorm;

                % recursively walk the channel path and unset anything
                % that is abandoned (i.e., cells that have no *other*
                % `flowsFrom` than the abandoned path.
                % note: this has to be recursive, because we may encounter branches
                % along the abandoned pathway and we need to clear each of
                % these pathways.
                abandonAll = false; % whether another channel stops abandonment or not
                [grid] = unmarkChannelToNode(grid, startIndex, prevIndex, abandonAll);

                % we removed one pathway (the jth pathway) of this branch,
                % but this changed the proportions and indices of channels
                % in the flowsToFracs array. We could devise some way to
                % reset this counter and keep going, but this would be
                % difficult, I think. Instead, we just stop and will catch
                % the abandonment on the next iteration.
                %break % break the `for j=` loop and go to next `i` branch

            end
        end
    end
end
