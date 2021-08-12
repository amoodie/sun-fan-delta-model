function grid = enactAvulsions(newAvulsions,grid,inlet)
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

    for n=1:numel(newAvulsions.rNew)

        % update distributary/tributary data in network
        % geometry
        grid.flowsTo{newAvulsions.indSource(n)} = [grid.flowsTo{newAvulsions.indSource(n)}; newAvulsions.indNew(n)]; % append new cell to the "flows to" list
        grid.flowsFrom{newAvulsions.indNew(n)} = [grid.flowsFrom{newAvulsions.indNew(n)}; newAvulsions.indSource(n)]; % append source cell to "flows from" list

        % start the new segment at the previously determined avulsion destination
        iStart = newAvulsions.rNew(n);
        jStart = newAvulsions.cNew(n);
        grid.channelFlag(iStart,jStart)=true; % flags new avulsion cell as channel
        % Propagate the avulsion channel if the current cell is
        % subaerial (i.e., oceanFlag == false)

        % check that all cells flagged as channels have defined
        % flowsFrom cells
        for i=1:grid.size(1)
            for j=1:grid.size(2)
                if  grid.channelFlag(i,j) && isempty(grid.flowsFrom{i,j}) && ne(i,inlet.row) && ne(j,inlet.col)
                    error('grid.flowsFrom not defined for channel cell');
                end
            end
        end

        if ~grid.oceanFlag(iStart,jStart)
            grid=propagateAvulsion(grid,iStart,jStart);
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