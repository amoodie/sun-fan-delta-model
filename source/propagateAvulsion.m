    function grid=propagateAvulsion(grid,avulsionCellInd) % nested function
       % propagateAvulsion: creates path for avulsion channel.
        %%% Deviation from the Sun et al. (2002) model:
        %%% their approach is weighted toward a steepest
        %%% descent, with randomness; however, the
        %%% implementation is ambiguous (particularly how
        %%% it affects path direction downstream of the
        %%% avulsion site). Instead implemented the
        %%% Karssenburg and Bridge (2008) approach that
        %%% just uses the steepest descent path (their p. 6); I start
        %%% the path from the detected avulsion destination
        %%% point. 
        %%%
        %%% To prevent single-cell loop and cross-over channel formation
        %%% during avulsion path finding, we restrict flow to certain
        %%% neighbor cells with the `forbiddenCells` array.

        % the current index is the avulsion cell
        indCurrent = avulsionCellInd;
        
        % The discharge at the branching cell.
        % After choosing the first step we will compute the approx
        % discharge partitioning at this new branch, and then use the
        % discharge down the new path to determine a critical stress for
        % stopping routing flow.

        % `indPrev` is used in the loop to ensure that flow never goes
        % directly back to where it came from. However, to initialize the
        % routing, we set the previous index as where the flow is already
        % going, so that this location cannot be the location of the first
        % step
        indPrev = grid.flowsTo{avulsionCellInd}';

        % configure index stepper based on grid dimensions
        iwalk = [-grid.size(1)-1, -1, +grid.size(1)-1, ...
                 +grid.size(1), +grid.size(1)+1, +1, -grid.size(1)+1, -grid.size(1)];

        % while there is still non-ocean non-channel non-sink cells to walk
        continuePropagateAvulsion = true;
        while continuePropagateAvulsion
            %% while we want to find the next step, find where *might* go

            % update list of forbiddenCells
            forbiddenCells = [avulsionCellInd, indPrev]; % previous location is forbidden

            % find the indices of the neighbors and get slopes to there
            nghbrs = indCurrent + iwalk;
            nghbrSlopes = [grid.S.NW(indCurrent) grid.S.N(indCurrent) grid.S.NE(indCurrent) ...
                           grid.S.E(indCurrent) grid.S.SE(indCurrent) grid.S.S(indCurrent) ...
                           grid.S.SW(indCurrent) grid.S.W(indCurrent)];
 
            % add to the list of forbidden cells with the locations that would create crossover channels
            [forbiddenCorners] = checkNeighborsChannelsCrossover(grid, nghbrs);
            forbiddenCells = unique([forbiddenCells, forbiddenCorners]);

            % prevent the flow from going to any cell that it already
            % receives flow from
            forbiddenCells = unique([forbiddenCells, grid.flowsFrom{indCurrent}']);

            % prevent the flow from going to any cell that is uphill
            nghbrSlopes(nghbrSlopes < 0) = NaN;

            % adjust slope array for the forbiddenCells, changes to NaN
            matches = ismember(nghbrs, forbiddenCells);
            nghbrSlopes(matches) = NaN;

            % do a safety check that some cell is finite
            % (choosable).
            if ~any(isfinite(nghbrSlopes))
                % break the while loop here instead of jsut setting route
                % to false, so that the *step is not taken*, i.e., no
                % avulsion is possible. In this scenario, the avulsion just
                % doesn't happen and the run can continue. Eventually an
                % avulsion will happen somewhere else and abandon this
                % pathway.
                break
            end

            % now find the index of the next location that we might visit
            [~,indSmax] = max(nghbrSlopes);

            %% take the step to determine what the new ind will be
            step = iwalk(indSmax);
            indNew = indCurrent + step;
            [iNew, jNew] = ind2sub(grid.size, indNew);

            % make the connection to this new cell
            grid.flowsTo{indCurrent} = [grid.flowsTo{indCurrent}; indNew]; % append new cell to the "flows to" list
            grid.flowsFrom{indNew} = [grid.flowsFrom{indNew}; indCurrent]; % append source cell to "flows from" list
            wasChannel = grid.channelFlag(indNew); % was this cell a channel *before* we got here
            grid.channelFlag(indNew) = true; % mark this cell as now being a channel

            %% determine whether the new point is somewhere we want to continue from

            % Stop path construction if new point is beyond domain boundary (Alternatively, could
            % have sidewalls steer flow (closed boundary) or make boundary open or periodic).
            if iNew<1 || iNew>grid.size(1) || jNew<1 || jNew==grid.size(2)

                continuePropagateAvulsion = false;

            % Stop path construction if new point was already marked as a
            % channel. In this scenario, we don't want to do any avulsion path
            % finding, the path is already determined.
            elseif wasChannel

                % end iteration of the while loop
                continuePropagateAvulsion = false;

            % Stop path construction if new point is an ocean cell. In this
            % scenario, we have reached an outlet, so stop routing.
            elseif grid.oceanFlag(indNew)

                % end iteration of the while loop
                continuePropagateAvulsion = false;

            % Stop path construction if new point is a sink cell. In this
            % scenario, we would only be able to go upslope to route of the
            % sink, and we want to stop pathfinding
            elseif grid.sinkFlag(indNew)

                % end iteration of the while loop
                continuePropagateAvulsion = false;

            % Continue path construction is the new point is an ordinary
            % land cell (non-ocean and non-channel cell)
            else
                % this is a land cell, so we need to find
                % a pathways across the land

                % make the currentInd this new step, to start over the
                % while loop and search for the next step to take
                indPrev = indCurrent;
                indCurrent = indNew;
            end

        end % end while loop for avulsion path construction 
    end % end nested function propagateAvulsion        


function [forbiddenCorners] = checkNeighborsChannelsCrossover(grid, nghbrs)
    %checkNeighborsChannelsCrossover Prevent crossover channels.
    %
    % Prevent the formation of crossover channels by checking neighboring
    % cells for channel connections and if they exist, using it to forbid
    % channels looping
    %
    % Examples of crossover channel cases to consider. In these examples,
    % imagine the avulsion path is being selected for channel 2, that is,
    % channel 1 already exists, and the next location chosen is the *.
    % Background (non-occupied) cells are #.
    %
    %   Case A      Case B      Case C
    %   -------     -------     -------
    %    1 1 *       2 # #       1 1 *
    %    2 2 1       1 2 1       2 2 3
    %    # # #       # 1 *       # 3 #
    %
    % Case A and Case B are physically impossible. These, we prevent from
    % forming by checking neighbors in a loop. Case C is potentially
    % possible, but may be unrealistic, given the spatial constraint and
    % channel widths. We use a parameter option flag to determine whether to
    % disallow this case (`cornerOption=='present'` to prevent these
    % avulsions, and `cornerOption=='connected'` to allow these avulsions).

    % Note in all cases, there are only four possible locations for the
    % crossover to form (the corners of the 3x3 grid). The approach here is
    % to define the cell pairs that constrain those locations, and loop
    % through them, determining if any cells should be prohibited.

    % constraining cell pairs, starting from N-E pair
    pairs = [2, 4; 4, 6; 6, 8; 8, 2]; % N-E; E-S; S-W; W-N
    corners = [3; 5; 7; 1]; % NW; NE; SE; SW

    % option flag
    cornerOption = 'connected'; % 'connected' | 'present'

    % preallocate the forbidden corners, shape is for four possible corners
    forbiddenCorners = false(1,4);

    % loop through
    for i=1:4
        ithPair = pairs(i,:);
        % always check the presence of channels first (faster)
        if all(grid.channelFlag(nghbrs(ithPair)))
            % both constraining cells are channels
            if strcmp(cornerOption, 'connected')
                % check whether these two constraining cells are actually
                % connected or are two different channels
                if any(grid.flowsTo{nghbrs(ithPair(1))} == nghbrs(ithPair(2))) || any(grid.flowsTo{nghbrs(ithPair(2))} == nghbrs(ithPair(1)))
                    % the channel cells are connected

                    % add this corner cell to the list of forbiddenCorners
                    forbiddenCorners(i) = true;
                end

            elseif strcmp(cornerOption, 'present')
                % doesn't matter if these are connected

                % add this corner to the list of forbiddenCorners
                forbiddenCorners(i) = true;

            end
        end
    end
    
    % convert the corners boolen to cell indices
    forbiddenCorners = nghbrs(corners(forbiddenCorners));
end
