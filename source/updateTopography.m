function grid = updateTopography(grid,inlet,lambda,tStep_sec,Qs_inlet)    
% updateTopography.m: updates bed elevation along flow paths using Equation 12.
    grid.Qs_in = zeros(grid.size);

    % determine all the channel locations to loop through
    channelInds = find(grid.channelFlag);

    % loop through the channel inds
    for kk=1:numel(channelInds)

        % grab this channelInd and conver to row,col
        k = channelInds(kk);
        [i,j] = ind2sub(grid.size,k);

        % look up the sediment discharge into the grid cell
        % If this is the first segment, then use inlet
        % sediment discharge. Otherwise, check the sediment
        % discharge to this cell using grid.flowsFrom
        if i==inlet.row && j==inlet.col
            grid.Qs_in(i,j) = Qs_inlet;
        else
            % look up the contributing cells
            contributingCells = grid.nghbrs(grid.flowsFromGraph(:,k), k);

            if numel(contributingCells)==1
                % one input, so takes all of output from contributor
                grid.Qs_in(i,j) = grid.Qs_out(contributingCells);
            elseif numel(contributingCells)>1 % multiple contributing cells
                % Sediment routing at network junctions is not explictly
                % discussed in the paper. For simplicitly, we assume
                % sediment discharge is partitioned in the same way 
                % that water discharge is partitioned.
                %   Implementation below is a vectorized computation for finding
                %   contributors and their relative sediment contributions
                Qs_outContributors = grid.Qs_out(contributingCells);  % Qs_out from each contributor
                channelContributorsBool = (grid.nghbrs(:, contributingCells) == k);  % bool where contributing cells are connected to kth cell
                flowsToFracContributors = grid.flowsToFrac(:, contributingCells);  % fracs for contributing cells connected to kth cell
                
                % determine the input from each contributor
                Qs_in_fromContributors = flowsToFracContributors(channelContributorsBool) .* Qs_outContributors;  % Qs into kth cell from contributing cells
                
                % update the field with sum of contributors
                grid.Qs_in(k) = sum(Qs_in_fromContributors);

            else
                % throw error if this cell is flagged as a channel
                % but there are no contributingCells defined,
                % unless this is the inlet cell (which has no
                % contributing cells)
                if grid.channelFlag(i,j) && ne(i,inlet.row) && ne(j,inlet.col)
                    error('No contributing cells found where they were expected');
                end
            end
        end
    end
    
    % set Qs_in, Qs_out to zero for non-channel cells (values are currently NaN)
    grid.Qs_in(~grid.channelFlag)=0;
    grid.Qs_out(~grid.channelFlag)=0;

    % conservation of mass (equation 12)
    grid.deltaz = ((grid.Qs_in-grid.Qs_out)/((1-lambda)*grid.dx^2))*tStep_sec; % equivalent to  -((Qs_out-Qs_in)/((1-lambda)*grid.dx^2))*tStep; 

    % update topography for all grid cells
    grid.z = grid.z + grid.deltaz;
end
