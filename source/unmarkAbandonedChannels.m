function grid=unmarkAbandonedChannels(grid,Qw_threshold)
    
    % make a list of all the branches
    channelStartIndices = [];
    for k=1:numel(grid.channelFlag)
        if numel(grid.flowsTo{k}) == 2
            channelStartIndices = [channelStartIndices, k];
        end
    end
    
    % while there are any branch indices remaining
    for i=1:numel(channelStartIndices)
        
        branchIndex = channelStartIndices(i);
        
        % try both pathways against the threshold discharge
        for j=1:numel(grid.flowsToFrac_Qw_distributed{branchIndex})
            testDischarge = grid.flowsToFrac_Qw_distributed{branchIndex}(j);
            if testDischarge < Qw_threshold
                % this channel pathways needs to be disconnected and unset
                % along the channel course. Get index of the step down the
                % abandoned channel
                startIndexList = grid.flowsTo{branchIndex}(j); % new path start
                prevIndexList = branchIndex; % previous step of new path start
                
                % unset flowsTo at the branch and update partitioning
                grid.flowsTo{branchIndex} = grid.flowsTo{branchIndex}(1:2 ~= j);
                grid.flowsToFrac_Qw_distributed{branchIndex} = 1;
                
                % now recursively walk the channel path and unset anything
                % that does not have a flowsFrom
                while numel(startIndexList) > 0
                    
                    startIndex = startIndexList(1);
                    prevIndex = prevIndexList(1);
                    startIndexList = startIndexList(2:end);
                    prevIndexList = prevIndexList(2:end);
                    
                    % recursively walk channels
                    [grid, newStarts, newPrevs] = unmarkChannelToNode(grid, startIndex, prevIndex);
                    startIndexList = [startIndexList; newStarts];
                    prevIndexList = [prevIndexList; newPrevs];

                end
                
                % if we removed one branch, we don't want to check the
                % other (i.e., skip the other j)
                break
            end
        end

    end
        
end

function [grid, newStarts, newPrevs] = unmarkChannelToNode(grid, startIndex, prevIndex)

    % set returns
    newStarts = []; % where to go
    newPrevs = []; % where came from

    % now walk the channel
    walk = true;
    currentIndex = startIndex;
    while walk

        if numel(grid.flowsFrom{currentIndex}) >= 2
            % this is another channnel, flowing from two places. Thus, it
            % will still have discharge, even after we abandon it from this
            % pathway. So, we want to remove the abandoned flowsFrom index,
            % and leave the loop without adding anything else to walk
            keepidx = (grid.flowsFrom{currentIndex} ~= prevIndex);
            grid.flowsFrom{currentIndex} = grid.flowsFrom{currentIndex}(keepidx);
            walk = false; % no more walking
            
        elseif numel(grid.flowsFrom{currentIndex}) <= 1
            % this is a single channel, so remove flowsFrom completely
            % for where we just left from
            grid.flowsFrom{currentIndex} = [];
            
            if numel(grid.flowsTo{currentIndex}) == 2
                % this cell is another branch
                grid.channelFlag(currentIndex) = false;
                
                % we need to walk both pathways, so return to the parent
                newStarts = grid.flowsTo{currentIndex};
                newPrevs = [currentIndex; currentIndex];
                
                % clear info on where this node would go
                grid.flowsTo{currentIndex} = [];
                
                walk = false; % no more walking here
            elseif numel(grid.flowsTo{currentIndex}) == 1
                % this is a simple channel
                % grab the next step
                nextIndex = grid.flowsTo{currentIndex};
                
                % now unset the flowsTo, channelFlag, partition
                grid.flowsTo{currentIndex} = [];
                grid.channelFlag(currentIndex) = false;
                grid.flowsToFrac_Qw_distributed{currentIndex} = [];
                
                prevIndex = currentIndex;
                currentIndex = nextIndex;
            elseif numel(grid.flowsTo{currentIndex}) == 0
                % this is an outlet
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
