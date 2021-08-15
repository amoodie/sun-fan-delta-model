    function grid=propagateAvulsion(grid,iStart,jStart,forbiddenCells) % nested function
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
        %%% To prevent single-cell loop formation during avulsion path
        %%% selection, we use the input array `forbiddenCells` to ensure that
        %%% the maximum slope direction is never the cell that flow is
        %%% avulsed from. To do this, we set any cell indices in
        %%% `forbiddenCells` equal to `NaN`, so that this index is not
        %%% selected during path finding.
       
        % extract the initial point into shorthand i, j
        i = iStart;
        j = jStart;
        indCurrent = sub2ind(grid.size, i, j);
        
        % configure index stepper based on grid dimensions
        iwalk = [-grid.size(1)-1, -1, +grid.size(1)-1, ...
                 +grid.size(1), +grid.size(1)+1, +1, -grid.size(1)+1, -grid.size(1)];

        % while there is still non-ocean non-channel cells to walk
        continuePropagateAvulsion = true;
        while continuePropagateAvulsion
            % find the indices of the neighbors and get slopes to there
            nghbrs = indCurrent + iwalk;
            nghbrSlopes = [grid.S.NW(indCurrent) grid.S.N(indCurrent) grid.S.NE(indCurrent) ...
                           grid.S.E(indCurrent) grid.S.SE(indCurrent) grid.S.S(indCurrent) ...
                           grid.S.SW(indCurrent) grid.S.W(indCurrent)];
            
            % adjust slope array for the forbiddenCells, changes to NaN
            matches = ismember(nghbrs, forbiddenCells);
            nghbrSlopes(matches) = NaN;
            
            % do a safety check that some cell is finite
            % (choosable).
            if ~any(isfinite(nghbrSlopes))
                error('all non-finite slopes encountered')
            end
            
            % now find the index of the next location to visit
            [Smax,indSmax] = max(nghbrSlopes);
    
            % now take the step
            step = iwalk(indSmax);
            indNew = indCurrent + step;
            [iNew, jNew] = ind2sub(grid.size, indNew);

           % Stop path construction if new point is beyond domain boundary (Alternatively, could
           % have sidewalls steer flow (closed boundary) or make boundary open or periodic). 
           
           if iNew<1 || iNew>grid.size(1) || jNew<1 || jNew==grid.size(2)
                continuePropagateAvulsion = false;
           else
                % update grid.flowsTo, grid.flowsFrom, grid.channelFlag using
                % the information for the current and next cell in the avulsion path
                grid.flowsTo{indCurrent} = [grid.flowsTo{indCurrent}; indNew];
                grid.flowsFrom{indNew} = [grid.flowsFrom{indNew}; indCurrent];
                grid.channelFlag(indNew) = true;

                % % check if the new cell intersects an existing channel or
                % ocean cell or topographic sink; if it does,
                % update flag to stop path construction. 
                if grid.channelFlag(indNew) || grid.oceanFlag(indNew) || grid.sinkFlag(indNew) 
                   continuePropagateAvulsion = false;
                else
                    % step forward by setting (i,j) to (iNew,jNew)
                    indCurrent = indNew;
                end
           end           
        end % end while loop for avulsion path construction 
    end % end nested function propagateAvulsion
