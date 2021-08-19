function grid = enactAvulsions(avulsionCellInds,grid,inlet)
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
    for n=1:numel(avulsionCellInds)

        % propogate the avulsion from the cell ind
        grid=propagateAvulsion(grid,avulsionCellInds(n));

        % check that all cells flagged as channels have defined
        % flowsFrom cells
        for i=1:grid.size(1)
            for j=1:grid.size(2)
                if  grid.channelFlag(i,j) && isempty(grid.flowsFrom{i,j}) && ne(i,inlet.row) && ne(j,inlet.col)
                    error('grid.flowsFrom not defined for channel cell');
                end
            end
        end

    end

    % check that all cells flagged as channels have defined
    % flowsFrom cells
    for i=1:grid.size(1)
        for j=1:grid.size(2)
            if  grid.channelFlag(i,j) && isempty(grid.flowsFrom{i,j}) && ne(i,inlet.row) && ne(j,inlet.col)
                error('grid.flowsFrom not defined for channel cell');
            end
        end
    end
end  % end path construction for avulsions