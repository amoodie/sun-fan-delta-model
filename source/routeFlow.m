function grid=routeFlow(grid,inlet,Qw_inlet,gamma,Qw_mismatch_tolerance)
    % routeFlow: route flow through the network
      %%% Flow routing begins with the inlet cell, and it formulated as a
      %%% sequential unordered walk through the channel network. Flow is
      %%% stepped cell-by-cell, only moving to flowsTo cells when a given
      %%% location has received flow from all contributing cells. This
      %%% algorithm is an O(n) process, where n is the number of cells in
      %%% the network.
      %%%
      %%% Routing is organized around `Qw_toRoute`, a `grid.size` of NaN, which
      %%% records the location and volume of flow in the network. The cell
      %%% to route flow next is selected by finding a location where
      %%% `~isnan(Qw_toRoute)` and `readyToFlow`, where `readyToFlow` is
      %%% true when a cell has received flow from all upstream cells. A
      %%% cell is also limited to only flow once, by checking if it has
      %%% `alreadyFlowed`; this check prevents loops in the network from
      %%% recursion, but may lead to flow becoming trapped within the
      %%% network. In these cases, the network loop is trimmed out, so that
      %%% there is no loop on the next timestep.

        %% set up the routing
        % routeFlow: routes flow along the channel network, splitting flow at branches using eqn. 11.
        grid.Qw = zeros(grid.size); % initialize grid that stores total water discharge during flow routing
        grid.Qw_toRoute = nan(grid.size); % initialize grid that stores water discharge to route
        grid.flowsToFrac(:) = 0;  % rezero grid to track fractional distribution of flow at distributaries
        alreadyFlowed = zeros(grid.size); % initialize grid to track where flow already visited (looping check)
        
        % To begin the routing, the only cell with flow is the inlet
        % cell; specify that discharge in Qw_grid
        grid.Qw_toRoute(inlet.row,inlet.col) = Qw_inlet;
        
        % set summed discharge for the inlet cell as well
        grid.Qw(inlet.row,inlet.col) = Qw_inlet;
        
        % get array index for inlet cell, the first cell from which to
        % route flow
        cellIndFlowToRoute = sub2ind(grid.size,inlet.row,inlet.col);

        % identify network endpoints as cells that do not flow into any
        % other cells
        grid.networkEndpoints = (grid.flowsToCount == 0);
        
        % determine the count of flows from
        grid.flowsFromCount = squeeze(sum(grid.flowsFromGraph, 1));
        
        % set the tracker for places that have all of their inputs filled
        grid.flowedFromSources = zeros(grid.size);
        grid.flowedFromSources(inlet.row,inlet.col) = 1;
        
        %% run the routing
        % while there are any cells that have flow to route
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
                
            % find where this cell is connected to
            channelsBool = grid.flowsToGraph(:, cellIndFlowToRoute);  % true where flowsTo has channels
            cellsIndFlowsTo = grid.nghbrs(channelsBool, cellIndFlowToRoute);  % the cell indices of the next cells

            if isempty(cellsIndFlowsTo)
                % doesn't flow anywhere
                %   set the current cell to NaN and leave
                grid.Qw_toRoute(cellIndFlowToRoute)=nan;
            elseif numel(cellsIndFlowsTo)==1
                % flows to one location
                %   set the next cell to have the current cell's total flow (plus whatever it already had!)
                grid.Qw_toRoute(cellsIndFlowsTo) = sum([grid.Qw_toRoute(cellsIndFlowsTo),Qw_routed],'omitnan');
                grid.Qw(cellsIndFlowsTo) = grid.Qw(cellsIndFlowsTo) + Qw_routed;
                % update partitioning information for current cell
                grid.flowsToFrac(channelsBool, cellsIndFlowsTo) = 1;
                % update flow information for current cell to empty
                grid.Qw_toRoute(cellIndFlowToRoute) = nan;
                grid.flowedFromSources(cellsIndFlowsTo) = grid.flowedFromSources(cellsIndFlowsTo) + 1;
            elseif numel(cellsIndFlowsTo)>=2
                % flow has reached a distributary junction

                % look up the slope to each cell from precomputed
                branchSlopes = grid.S.d8(channelsBool, cellIndFlowToRoute);
                if all(branchSlopes < 0)
                    % all negative
                    branchFractions = (abs(branchSlopes).^-gamma) / sum(abs(branchSlopes).^-gamma);
                else
                    % may be all positive, or some mix of positive and negative and 0
                    %  Routing is according to Murray and Paola:
                    %  > The sum, which runs over all the neighbours with 
                    %  > positive slopes, normalizes the water routing
                    %  so we set any negative slopes to 0 so that they
                    %  receive no flow
                    branchSlopes(branchSlopes < 0) = 0; % set <0 to 0
                    branchFractions = (branchSlopes.^gamma) / sum(branchSlopes.^gamma);
                end

                % record the computed flow partitioning
                grid.flowsToFrac(channelsBool, cellIndFlowToRoute) = branchFractions;

                % check that the routed discharge is positive for both
                % branches, and that sum of discharge going to the two
                % distributaries equals the discharge from the contributing
                % cell
                if any(branchFractions < 0)
                    error('Negative discharge routed at distributary junction')
                elseif abs(grid.Qw_toRoute(cellIndFlowToRoute) - (Qw_routed * sum(branchFractions))) > (Qw_mismatch_tolerance)
                    error('Inconsistent mass balance for discharge at distributary junction')
                end

                % update grid.Qw and grid.Qw_toRoute
                grid.Qw(cellsIndFlowsTo) = grid.Qw(cellsIndFlowsTo) + (Qw_routed * branchFractions);
                grid.Qw_toRoute(cellsIndFlowsTo) = sum([grid.Qw_toRoute(cellsIndFlowsTo), Qw_routed * branchFractions],2,'omitnan');
                grid.Qw_toRoute(cellIndFlowToRoute) = NaN;
                
                % record the flow having gone into each of the cells
                grid.flowedFromSources(cellsIndFlowsTo) = grid.flowedFromSources(cellsIndFlowsTo) + 1;
                
            else
                error('Flow routing: Invalid value for cellsIndFlowsTo');
            end

            % update the status of the loop checking array for the current cell
            alreadyFlowed(cellIndFlowToRoute) = 1;

            % update list of cells from which to route flow for next iteration
            %   find where there is water to route
            hasFlow = ~isnan(grid.Qw_toRoute);
            % find where cells have received flow from all upstream cells
            hasSources = (grid.flowedFromSources == grid.flowsFromCount);
            % find where a cell is ready to flow
            readyToFlow = and(hasFlow,hasSources);
            % places that haven't already been flowed (prevents loops)
            notAlreadyFlowed = ~(alreadyFlowed);  % prevents loops

            % choose the next cell to flow as any one that is ready to
            %    flow, but hasn't already flowed
            cellIndFlowToRoute = find(and(readyToFlow,notAlreadyFlowed),1,'first');

        end

        %% flow routing is now complete. Run some diagnostic checks:
        % (1) check if there is any flow remaining in the domain
        if any(~isnan(grid.Qw_toRoute(:)))
            % if there is any flow remaining, then there is a loop.
            %    we trim out loops here by walking down stream from the
            %    loop until encountering another channel pathway.
            %    Note: we walk down all paths that the cell flowsTo,
            %    because we do not know the loop structure.

            % find remaining water
            remaining_Qw = find(~isnan(grid.Qw_toRoute));

            % recursively walk the channel path and unset channels
            % note: this function is recursive. See docstring and
            % `unmarkAbandonedChannels.m` for more information.
            abandonAll = false; % whether another channel stops abandonment or not
            % for each remaining cell
            for i = 1:length(remaining_Qw)
                remaining_Qw_i = remaining_Qw(i);
                % find where it flows to, then remove connection to
                channelsBool = grid.flowsToGraph(:, remaining_Qw_i);
                flowsTo_i = grid.nghbrs(channelsBool, remaining_Qw_i);
                grid.flowsToCount(remaining_Qw_i) = 0;
                grid.flowsToGraph(:, remaining_Qw_i) = 0;
                % for each place it flows to
                for j = 1:length(flowsTo_i)
                    % unset the channels down the path
                    [grid] = unmarkChannelToNode(grid, flowsTo_i(j), remaining_Qw(i), abandonAll);
                end
            end
        else
            % if there is no flow remaining
            % (1a) sum the discharge for the network endpoints and verify that
            %    it equals the input discharge
            Qw_out_total = sum(grid.Qw(grid.networkEndpoints));
            if abs(Qw_out_total - Qw_inlet) > Qw_mismatch_tolerance
                error('flowRoute: Summed discharge at channel network outlets does not equal discharge at inlet');
            end
        end

        % (2) Check that no cell has a discharge greater than the inlet discharges 
        if any((grid.Qw(:) - Qw_inlet) > Qw_mismatch_tolerance)
            error('Error in flow routing: discharge greater than input discharge detected')
        end

    end % end nested function routeFlow
