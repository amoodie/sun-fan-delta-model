function [grid] = validateNetwork(grid,inlet)
% validateNetwork Determine the network is a acyclic directed graph.
    % This function determines whether there are any *invalidating* loops,
    % and removes any that are found. The algorithm is recursive,
    % essentially walking down each flow pathway encountered, and keeping a
    % record of everywhere we have been (prevIndList).
    %
    % In the end, the grid network returned is a valid acyclic directed
    % graph, meaning that there is one source node (inlet) and all pathways
    % flow to an outlet without returning to where they have already been.
    % This allows for local loops, e.g., where flow branches and later
    % converges downstream and continues to flow to an outlet.

    % allocate the "where we have been list" as empty
    prevIndList = [];

    % start the walking from the inlet
    startInd = sub2ind(grid.size, inlet.row, inlet.col);
    
    % walk.
    % this algorithm is recursive, so the same function will be called when
    % any branches are encountered.
    doTrim = true; % whether to actually trim (or just search)
    [grid, ~] = walkToNodeSearchForLoops(grid, startInd, prevIndList, doTrim);
    
end

function [grid, trimFlag] = walkToNodeSearchForLoops(grid, startInd, prevIndList, doTrim)
% walkToNodeSearchForLoops walks down a pathway, finding and trimming loops
    % This function does the actual walk down the pathways (channels)
    % looking for loops, i.e., places where the pathway leads back to
    % somewhere it has already been. Note, however that local branching and
    % convergence is allowed.
    %
    % Implementation is recursive, so is called again when a branch is
    % encountered.
    %
    % Inputs:
    %   grid
    %   startInd - where the pathway starts
    %   prevIndList - all the places we have already been on this pathway
    %   doTrim - boolean for whether to actually trim, or just search for loops

    trimFlag = false; % whether a loop is encountered and trimming is started
    step = true; % whether to continue stepping down the path

    % update the list of previous places we have been with the starting
    % point of this path
    prevIndList = [prevIndList; startInd];

    % update our current location with the starting point of this path
    currentInd = startInd;

    while step

        % where the next step will be to
        nextInd = grid.flowsTo{currentInd};

        % if the next step is anywhere we have been
        if ismember(nextInd, prevIndList)
            % A loop has been found, so we have to trim
            trimFlag = true; % flag that a loop has been found
            if doTrim % if we actually want to do the trimming
                % call a subfunction with routines to trim out the loop
                [grid] = trimOutLoop(grid, nextInd, prevIndList);
            end
            % now break the while loop
            % *without proceeding to the next if statement*
            break
        end

        % if the next step were to be into a loop, we do not reach this
        % point. So, since there is not a loop at the next step, we need to
        % determine what type of cell the next step is, and then respond
        % accordingly.
        if numel(nextInd) == 2
            % the next step is a branching cell

            % update the list of places we have been with the current ind
            %%% I actually don't think this is necesary, because it should
            %%% have been taken care of in either stepping below or in
            %%% start of this function for sub-branches.
            prevIndList = [prevIndList; currentInd];
            
            % walk down each pathway
            for i=1:2
                % walk each path, recursively
                doTrim = true;
                [grid, trimFlag] = walkToNodeSearchForLoops(grid, nextInd(i), prevIndList, doTrim);
            end
            % now can stop this loop since all below walked
            step = false; 

        elseif numel(nextInd) == 1
            % the next step is a conduit channel.
            % take a step
            
            % update the list of where we have been with the current cell
            prevIndList = [prevIndList; currentInd];
            
            % take the step to the next cell and restart the while loop
            currentInd = nextInd;

        elseif numel(nextInd) == 0
            % the next step is an outlet, can stop
            step = false;

        else
            error('Error in number of indices of next step')

        end
    end
end

function [grid] = trimOutLoop(grid, nextInd, prevIndList)
    % we want to unset the channel cells between *here* and the
    % last place that flows *somewhere other than here*.
    % To find that place, we walk *up* the previous index list, and
    % look for the first branch that does not return a trim
    % flag (i.e., it flows somewhere else).
    trimmed = false;
    for k=0:numel(prevIndList)-1
        checkInd = prevIndList(end-k);
        if numel(grid.flowsTo{checkInd}) == 1
            % not a branch, no need to check
            continue
        elseif numel(grid.flowsTo{checkInd}) == 0
            % this should be impossible...
            disp('numel was 0, what does this mean? Exit loop not working?');
        elseif numel(grid.flowsTo{checkInd}) == 2
            % this is a branch, check the other pathway
            thisPathInd = prevIndList(end-k+1);
            otherPathInd = grid.flowsTo{checkInd}(grid.flowsTo{checkInd} ~= thisPathInd);
            doTrim = false;
            [~, trimFlag] = walkToNodeSearchForLoops(grid, otherPathInd, prevIndList, doTrim);
            if trimFlag
                % if this pathway still leads to a loop
                % continue walking up the path to the next branch
                continue
            else
                % this is a branch we can split the network at!
                grid.flowsTo{checkInd} = [otherPathInd];
                grid.flowsToFrac_Qw_distributed{checkInd} = 1;

                % we recursively walk down this path to unmark it
                startIndex = thisPathInd; % the head of the path to abandon
                prevIndex = checkInd; % the branching cell
                [grid] = unmarkChannelToNode(grid, startIndex, prevIndex);

                % now we can exit the for loop
                trimmed = true;
                break
            end
        end
    end

    % if trimmed is false, we did not find a thing to trim, need to debug
    if ~trimmed
        % no trimmable loop was found at a branch, so we just trim 
        % everything downstream of the looped cell index.

        % find the index of the nextInd in the prevIndexList
        whr = find(nextInd == prevIndList,1);

        % unset the flow to/from this location
        grid.flowsTo{nextInd} = [];
        grid.flowsToFrac_Qw_distributed{nextInd} = [];

        % find the pathway head
        head = min([whr+1; length(prevIndList)]);
        startIndex = prevIndList(head); % the head of the path to abandon
        prevIndex = nextInd; % the branching cell
        [grid] = unmarkChannelToNode(grid, startIndex, prevIndex);
    end
end