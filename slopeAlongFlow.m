
    function grid=slopeAlongFlow(grid)
        % slopeAlongFlow: retrieves cell-to-cell slope values along flow
        % paths
        
        grid.S.alongFlow = nan(grid.size); 
     
        for k=1:numel(grid.channelFlag)
            if grid.channelFlag(k) 
                if isempty(grid.flowsTo{k})
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
                    S_temp = nan(size(flowsTo));
                    for l=1:numel(flowsTo)
                        S_temp(l) = -(grid.z(flowsTo(l)) - grid.z(k)); % downhill slopes are defined as positive
                    end
                    grid.S.alongFlow(k) = max(S_temp); % if there are multiple flowsTo cells, treat slope as the maximum of slopes to each cell
                end         
            end
        end        
    end % end nested function slopeAlongFlow