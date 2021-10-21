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
            contributingCells = grid.flowsFrom{i,j};

            if numel(contributingCells)==1
                grid.Qs_in(i,j) = grid.Qs_out(contributingCells);
            elseif numel(contributingCells)>1 % multiple contributing cells
                for cc=1:numel(contributingCells)
                    % Sediment routing at network junctions is not explictly
                    % discussed in the paper. For simplicitly, I assume
                    % that the sediment discharge is partitioned in the
                    % same way that water discharge is partitioned.

                    % look up Qs_out at the contributing cell; use the water discharge partition fraction at
                    % at the tributary to set the partition fraction
                    % for sediment
                    Qs_out_fromContributingCell = grid.Qs_out(contributingCells(cc));
                    % Qs_out_fromContributingCell is sometimes empty:
                    % if discharge did not reach this location due to routeFlow rules,
                    % depsite these cells being connected by a channel in the tracking cell arrays.
                    % We thus include a safety here, to just proceed with the Qs calculation,
                    % but with 0 sediment contribution from this cell
                    if Qs_out_fromContributingCell > 0
                        Qw_fraction_receivedFromContributingCell = grid.flowsToFrac_Qw_distributed{contributingCells(cc)};
                        [receivingCellRow,receivingCellCol] = ind2sub(grid.size,grid.flowsTo{contributingCells(cc)});
                        if i==receivingCellRow(1) && j==receivingCellCol(1)
                            Qw_fraction_receivedFromContributingCell = Qw_fraction_receivedFromContributingCell(1);
                        elseif i==receivingCellRow(2) && j==receivingCellCol(2)
                            Qw_fraction_receivedFromContributingCell = Qw_fraction_receivedFromContributingCell(2);
                        else
                            error('Error in data lookup for contributing cell')
                        end
                        grid.Qs_in(i,j) = grid.Qs_in(i,j) + Qs_out_fromContributingCell*Qw_fraction_receivedFromContributingCell;
                    end
                end
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

    % set Qs_in, Qs_out to zero for non-channel cells (values are
    % currently NaN)
    grid.Qs_in(~grid.channelFlag)=0;
    grid.Qs_out(~grid.channelFlag)=0;

    % conservation of mass (equation 12)
    grid.deltaz = ((grid.Qs_in-grid.Qs_out)/((1-lambda)*grid.dx^2))*tStep_sec; % equivalent to  -((Qs_out-Qs_in)/((1-lambda)*grid.dx^2))*tStep; 

    % update topography for all grid cells
    grid.z = grid.z + grid.deltaz;
end
