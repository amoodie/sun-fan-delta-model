function grid=validateNetwork(grid,inlet)

    prevIndList = [];

    startInd = sub2ind(grid.size, inlet.row, inlet.col);
    
    % walk recursively
    [grid] = walkToNode(grid, startInd, prevIndList, 0);
    
end

function [grid, trimFlag] = walkToNode(grid, startInd, prevIndList, lvl)

    trimFlag = false;

    step = true;
    prevIndList = [prevIndList; startInd];
    currentInd = startInd;
    while step
        
        % next step
        nextInd = grid.flowsTo{currentInd};
        
        % if the next step is anywhere we have been
        if ismember(nextInd, prevIndList)
            % we have to trim
            trimFlag = true;
            
            % we want to unset the channel cells between *here* and the
            % last place that flows *somewhere other than here*.
            % To find that place, we walk *up( the previous index list, and
            % look for the first branch that does not return a trim
            % flag (i.e., it flows somewhere else).
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
                    [~, trimFlag] = walkToNode(grid, otherPathInd, prevIndList, lvl+1);
                    if trimFlag
                        % if this pathway still leads to a loop
                        % continue walking up the path to the next branch
                        continue
                    else
                        % this is a branch we can split the network at!
                        grid.flowsTo{checkInd} = [otherPathInd];
                        grid.flowsToFrac_Qw_distributed{checkInd} = 1;
                        
                        %% NEED TO ADD THE CODE TO ABANDON HERE!!!
                        % we recursively walk down this path to unmark it
                        startIndex = thisPathInd; % the head of the path to abandon
                        prevIndex = checkInd; % the branching cell
                        [grid] = unmarkChannelToNode(grid, startIndex, prevIndex);
                        
                        % now we can exit the for loop
                        break
                        % should we be leaving the entire walk? Do I want
                        % to be doing this abandoning in another func? Why
                        % do we end up in the numel==0 pathway after this?
                    end
                end
            end
        end

        % determine the next step type
        if numel(nextInd) == 2
            % it is a branch.
            prevIndList = [prevIndList; currentInd];
            for i=1:2
                % walk each path, recursively
                [grid, trimFlag] = walkToNode(grid, nextInd(i), prevIndList, lvl+1);
            end
            step = false; % now can stop this loop since all below walked

        elseif numel(nextInd) == 1
            % it is a conduit channel.
            % take a step
            prevIndList = [prevIndList; currentInd];
            currentInd = nextInd;

        elseif numel(nextInd) == 0
            % outlet, can stop
            step = false;

        else
            error('Error in number of indices of next step')

        end
    end
end