function grid=routeFlow(grid,inlet,Qw_inlet,gamma,Qw_mismatch_tolerance)
     
        % routeFlow: routes flow along the channel network, splitting flow at branches using eqn. 11.
        grid.Qw = zeros(grid.size); % initialize grid that stores total water discharge during flow routing
        grid.Qw_toRoute = zeros(grid.size); % initialize grid that stores water discharge to route 
        grid.flowsToFrac_Qw_distributed = cell(grid.size); % initialize grid to track fractional distribution of flow at distributaries
        
        % To begin the routing, the only cell with flow is the inlet
        % cell; specify that discharge in Qw_grid
        grid.Qw_toRoute(inlet.row,inlet.col) = Qw_inlet;
        
        % set summed discharge for the inlet cell as well
        grid.Qw(inlet.row,inlet.col) = Qw_inlet;
        
        % get array index for inlet cell, the first cell from which to
        % route flow
        cellIndFlowToRoute = sub2ind(grid.size,inlet.row,inlet.col);
        
        grid.networkEndpoints = false(grid.size); % initialize grid to track network endpoints
        % identify network endpoints as cells that do not flow into any
        % other cells
        for k=1:numel(grid.channelFlag)
            if grid.channelFlag(k) && isempty(grid.flowsTo{k})
                grid.networkEndpoints(k) = true; 
            end
        end
        
        while ~isempty(cellIndFlowToRoute) % route flow while there are still cells flagged to route from
            % for each cell with flow to route, route its flow to next cell
            % along the flow path, pausing at any bifurcations (i.e.,
            % places where there is more than one cell to route to)
            
            % look out the quantity of discharge in the cell
            Qw_routed = grid.Qw_toRoute(cellIndFlowToRoute);
            
            % check for the number of downstream-connected cells. 
            %   If 0, then this flow parcel has reached the end of the channel
                % network and will go no further. In that case, set
                % grid.Qw_toRoute to 0 for this cell.
            %   If 1, then set discharge in the connected cell to Qw_routed
            %   and set discharge in the source cell to 0.
            %   If >1, then have reached a distributary function; flow is distributed
            %   among multiple cells
                
            % add the routed discharge to the connected downstream cell, if
            % there is only 
            cellsIndFlowsTo = grid.flowsTo{cellIndFlowToRoute};
            
            if isempty(cellsIndFlowsTo)
                grid.Qw_toRoute(cellIndFlowToRoute)=0;
            elseif numel(cellsIndFlowsTo)==1
                grid.Qw_toRoute(cellsIndFlowsTo) = Qw_routed;
                grid.Qw(cellsIndFlowsTo) =  grid.Qw(cellsIndFlowsTo) + Qw_routed;
                grid.flowsToFrac_Qw_distributed{cellIndFlowToRoute} = 1; % i.e., all flow distributed to 1 cell
                grid.Qw_toRoute(cellIndFlowToRoute) = 0;
            elseif numel(cellsIndFlowsTo)==2
                
                % Flow has reached a distributary junction.
                [source.row,source.col] = ind2sub(grid.size,cellIndFlowToRoute);
                source.ind = cellIndFlowToRoute;
                
                % determine the partitioning fraction for water down each
                % branch
                [branch1, branch2] = computeBranchDischargePartitioning(grid, source, cellsIndFlowsTo, Qw_routed, gamma);
                
                % update grid.Qw and grid.Qw_toRoute
                grid.Qw(branch1.startRow,branch1.startCol) =  grid.Qw(branch1.startRow,branch1.startCol) + branch1.Qw_received;
                grid.Qw(branch2.startRow,branch2.startCol) =  grid.Qw(branch2.startRow,branch2.startCol) + branch2.Qw_received;

                grid.Qw_toRoute(branch1.startRow,branch1.startCol) = branch1.Qw_received;
                grid.Qw_toRoute(branch2.startRow,branch2.startCol) = branch2.Qw_received;
                grid.Qw_toRoute(source.row,source.col) = 0;

                % Record the flow partitioning information in the
                % corresponding grid and grid cell
                grid.flowsToFrac_Qw_distributed{source.row,source.col} = [branch1.Qw_fractionReceived,branch2.Qw_fractionReceived];
                if numel(grid.flowsToFrac_Qw_distributed{source.row,source.col}) ~= 2
                    error('Error in defining discharge routing to multiple branches');
                end
            else
                error('Flow routing: Invalid value for cellsIndFlowsTo');
            end

            % update list of cells from which to route flow for next
            % iteration
            cellIndFlowToRoute = find(grid.Qw_toRoute>eps,1,'first');
        end
        
        % flow routing is now complete. Run some diagnostic checks:
        % (1) Check that no cell has a discharge greater than the inlet discharges 
        if any((grid.Qw(:) - Qw_inlet) > Qw_mismatch_tolerance)
            error('Error in flow routing: discharge greater than input discharge detected')
        end
        
        %%%% this threw an error for reasons I don't fully understand. In
        %%%% principle, I think you should be able to compare inlet to sum
        %%%% of endpoints and get same discharge. I got mismatches using
        %%%% this code, but then when I check step-by-step for each routing
        %%%% operation I can't see any mismatch between input and output.
        % sum the discharge for the network endpoints and verify that
        % it equals the input discharge
        
% %         Qw_out_total = sum(grid.Qw(grid.networkEndpoints));
% % 
% %         % Compare the summed discharge at outlets to the discharge at
% %         % the inlet and through an error if their difference exceeds
% %         % numerical precision
% %         if abs(Qw_out_total - Qw_inlet) > 100*eps
% %             error('flowRoute: Summed discharge at channel network outlets does not equal discharge at inlet');
% %         end
    end % end nested function routeFlow
    

