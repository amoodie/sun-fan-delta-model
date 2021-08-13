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

        % first we want to check if this cell is already marked as a
        % channel; in this scenario, we don't want to do any avulsion path
        % finding, the path is already determined. This scenarios arises
        % when 1) we avulse into a channel that simply exists in a
        % neighboring cell, but it is also arises when multiple avulsions
        % are identified in a single `avulsionCheck` call, several of which
        % might flow into the same cell.
        if grid.channelFlag(iStart, jStart)
            % if this cell is already a channel, we don't want to do any
            % avulsion path finding.
            continue % end this iteration of the for loop
        elseif grid.oceanFlag(iStart,jStart)
            % this is an ocean cell, we don't want to do any avulsion
            % path finding
            grid.channelFlag(iStart,jStart)=true; % flags new avulsion cell as channel
            continue % end this iteration of the for loop
        else
            % this is a non-ocean and non-channel cell, so we need to find
            % a pathways across the land
            grid.channelFlag(iStart,jStart)=true; % flags new avulsion cell as channel

            % start the path finding for the new avulsion-created channel.
            % This process starts from the new avulsion-created channel head
            % (iStart, jStart), and not from the avulsion site.
            grid=propagateAvulsion(grid,iStart,jStart,[newAvulsions.indSource]);
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