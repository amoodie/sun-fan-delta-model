function [grid,successfulAvulsions] = enactAvulsions(avulsionCellInds,grid,inlet,branchLimit,propogationRule)
% enactAvulsions.m: Generates new channel flow paths based on new avulsion
% sites identified with avulsionCheck.m. The new channel path originating 
% from each avulsion site is constructed following paragraph 20 and equation 15.
% The path construction stops when any of the following conditions are met:
% (1) path meets an existing channel, (2) path reaches standing water (a 
% "wet" cell), or (3) path runs into a topographic sink (local minimum). 

    % identify sinks as cells with a difference between the filled and
    % unfilled DEM greater than numerical precision. Sinks are used to
    % halt the consturction of new avulsion paths. 
    grid.zFill = imfill(grid.z);
    grid.sinkFlag = (grid.zFill - grid.z) > eps;

    % loop through each cell identified for avulsion
    numAvulsions = size(avulsionCellInds, 1);
    successfulAvulsions = NaN(size(avulsionCellInds,1));
    for navul=1:numAvulsions

        % navul avulsion location
        navulInd = avulsionCellInds(navul,:);

        if ~grid.channelFlag(navulInd(1))
            break
        end

        % propogate the avulsion from the cell ind
        [grid,navulFirstStep] = propagateAvulsion(grid,navulInd,propogationRule);

        % if the avulsion was successful, add it to the list of avulsions
        % to return to be recorded
        if (navulFirstStep>0)
            successfulAvulsions(navul) = navulInd(1);
        end

        % if the number of branches is above the branch threshold, we have to
        % abandon one of the old channels.
        if (navulFirstStep>0) && (branchLimit < sum(grid.flowsToGraph(:, navulInd(1))))

            % choose one of the old branches to unmark
            allBranches = grid.nghbrs(grid.flowsToGraph(:, navulInd(1)), navulInd(1));
            oldBranches = allBranches(allBranches~=navulFirstStep);
            % for now, choose one at random
            startIndex = oldBranches(randi([1, branchLimit], 1));
            prevIndex = navulInd(1);  % the branch location
            branchIndex = navulInd(1);

            % unset flowsTo *at the branch* and update partitioning
            which_index = find((startIndex - branchIndex) == grid.iwalk);
            grid.flowsToGraph(which_index, branchIndex) = 0;
            grid.flowsToFrac(which_index, branchIndex) = 0;

            % renormalize the other outputs to equal 1
            renorm = grid.flowsToFrac(:, branchIndex) / sum(grid.flowsToFrac(:, branchIndex));
            grid.flowsToFrac(:, branchIndex) = renorm;

            abandonAll = false; % whether another channel stops abandonment or not
            [grid] = unmarkChannelToNode(grid, startIndex, prevIndex, abandonAll);
        end
        % determine the count of flows from
        grid.flowsFromCount = squeeze(sum(grid.flowsFromGraph, 1));
        grid = countFlowsToInds(grid);
        
    end

    % clean up the successful avulsions record
    successfulAvulsions = successfulAvulsions(~isnan(successfulAvulsions));

    % check that all cells flagged as channels have defined
    %    flowsFrom cells
    channelInds = find(grid.channelFlag);
    % unmark the inlet
    inlet_ind = sub2ind(grid.size, inlet.row, inlet.col);
    channelInds(channelInds == inlet_ind) = [];
    for kk=1:numel(channelInds)
        k = channelInds(kk);
        if  grid.channelFlag(k) && (sum(grid.flowsFromGraph(:,k)) == 0)
            error('grid.flowsFrom not defined for channel cell');
        end
    end

end  % end path construction for avulsions
