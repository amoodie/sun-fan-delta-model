
    function grid=propagateAvulsion(grid,iStart,jStart) % nested function
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
       
        i = iStart;
        j = jStart;
        
        continuePropagateAvulsion = true;
        while continuePropagateAvulsion
            % identify neighboring cell associated with maximum slope
            [Smax,indSmax] = max([grid.S.NW(i,j) grid.S.N(i,j) grid.S.NE(i,j) ...
               grid.S.E(i,j) grid.S.SE(i,j) grid.S.S(i,j) grid.S.SW(i,j) grid.S.W(i,j)]);

            if Smax == 0 % this if condition is no longer needed now that have imparted closed boundaries
               iNew = i+1; % extend channel toward bottom of grid (advance by 1 row)
               jNew = j;
            else
               if indSmax == 1 % NW
                   iNew = i-1;
                   jNew = j-1;
               elseif indSmax == 2 % N
                   iNew = i-1;
                   jNew = j;
               elseif indSmax == 3 % NE
                   iNew = i-1;
                   jNew = j+1;
               elseif indSmax == 4 % E
                   iNew = i;
                   jNew = j+1;
               elseif indSmax == 5 % SE
                   iNew = i+1;
                   jNew = j+1;
               elseif indSmax == 6 % S
                   iNew = i+1;
                   jNew = j;
               elseif indSmax == 7 % SW
                   iNew = i+1;
                   jNew = j-1;
               elseif indSmax == 8 % W
                   iNew = i;
                   jNew = j-1;
               end
            end

           % Stop path construction if new point is beyond domain boundary (Alternatively, could
           % have sidewalls steer flow (closed boundary) or make boundary open or periodic). 
           
           if iNew<1 || iNew>grid.size(1) || jNew<1 || jNew==grid.size(2)
                continuePropagateAvulsion = false;
           else
                % update grid.flowsTo, grid.flowsFrom, grid.channelFlag using
                % the information for the current and next cell in the avulsion path
                indCurrent = sub2ind(grid.size,i,j);
                indNew = sub2ind(grid.size,iNew,jNew);                    
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
                    i = iNew;
                    j = jNew;
                end
           end           
        end % end while loop for avulsion path construction 
    end % end nested function propagateAvulsion        