    function grid=slopeAlongFlow(grid)
        % slopeAlongFlow: retrieves cell-to-cell slope values along flow
        % paths

        % preallocate with NaNs
        grid.S.alongFlow = nan(grid.size); 

        % determine all the channel locations to loop through
        channelInds = find(grid.channelFlag);

        for kk=1:numel(channelInds)
            k = channelInds(kk);
            if grid.flowsToCount(k) == 0
                grid.S.alongFlow(k) = 0; % define slope as zero for any cell that does not flow elsewhere. 
                % Other conditions should kick into effect to treat
                % these cells, namely: (1) Qs_out is set to zero for all
                % cells that don't flow to other cells; and (2) a zero
                % value of slope causes an undefined value of flow
                % depth, and flow depth influences the avulsion
                % criterion (eqn. 13). For cases with negative/NaN flow
                % depths, the depth used in the avulsion criterion
                % (eqn. 13) is set to zero so that it does not affect
                % avulsion susceptibility.
            else
                % if this cell is a channel cell, look up the cell(s) it
                % flows to
                flowsTo = grid.flowsTo{k};
                [ri,ci] = ind2sub(grid.size, k); % row,col of cell i
                S_temp = nan(size(flowsTo));
                for l=1:numel(flowsTo)
                    [rj,cj] = ind2sub(grid.size, flowsTo(l)); % row,col of cell j
                    distance = sqrt((rj-ri)^2 + (cj-ci)^2) * grid.dx;  % distance along flow
                    S_temp(l) = -(grid.z(flowsTo(l)) - grid.z(k)) / distance; % downhill slopes are defined as positive
                end
                grid.S.alongFlow(k) = max(S_temp); % if there are multiple flowsTo cells, treat slope as the maximum of slopes to each cell
            end
        end
    end % end nested function slopeAlongFlow
