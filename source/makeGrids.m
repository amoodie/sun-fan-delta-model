function grid =  makeGrids(grid,oceanLevel) % nested function
    % makeGrids: generates elevation and other grids.   
     
        grid.xVec = 0:grid.dx:grid.xExtent; % vector of x-coordinates of grid, m
        grid.yVec = (grid.yExtent:-grid.dx:0)';% vector y-coordiantes of grid, m
        [grid.x,grid.y] = meshgrid(grid.xVec,grid.yVec); % makes 2D arrays of x- and y- coordinates of grid, m
        grid.size = size(grid.x);
        
        switch grid.DEMoptions.initialSurfaceGeometry.type
            case 'flat' % flat surface with no roughness
                grid.z = grid.DEMoptions.initialSurfaceGeometry.minElev*ones(grid.size);
            case 'slopeBreak'
                % Create two sections of equal length, each with constant slope
                rowSlopeBreak = round(grid.size(1)*grid.DEMoptions.slopeBreak.fracDistToSlopeBreak);
                slopeVec = zeros(grid.size(1),1);
                slopeVec(1:rowSlopeBreak) = grid.DEMoptions.slopeBreak.slope1;
                slopeVec(rowSlopeBreak+1:end) = grid.DEMoptions.slopeBreak.slope2;
                zVec(1)=0;
                for i=2:grid.size(1)
                    zVec(i,1)=zVec(i-1,1)+grid.dx*slopeVec(i);
                end
                grid.z = repmat(zVec,1,grid.size(2));
                grid.z = grid.z - (min(grid.z(:))-grid.DEMoptions.initialSurfaceGeometry.minElev);

                % Record row of slope break and corresponding elevation
                grid.DEMoptions.rowSlopeBreak = rowSlopeBreak;
                grid.DEMoptions.elevSlopeBreak = zVec(rowSlopeBreak);
                
                if grid.DEMoptions.slopeBreak.carveChannel
                    rows = 1:rowSlopeBreak; % carves channel on upstream side of each slope break in periodic topography
                    col = round(grid.size(2)/2); % middle column of DEM
                    channelDepthColumn = zeros(numel(rows),1);
                    channelDepthColumn(1:rowSlopeBreak) = grid.DEMoptions.slopeBreak.channelDepth;
                    % taper smoothly to zero
                    len = round(grid.size(1)*grid.DEMoptions.slopeBreak.channelDepthTaperLength);
                    channelDepthTaperDownstream = linspace(grid.DEMoptions.slopeBreak.channelDepth,0,len)';
                    channelDepthColumn((rowSlopeBreak-len+1):rowSlopeBreak) = channelDepthTaperDownstream;
                    grid.z(rows,col)= grid.z(rows,col) - channelDepthColumn;
                end
        end
        
        if grid.DEMoptions.addNoise
            % Add random noise 
            grid.z = grid.z + grid.DEMoptions.noiseAmplitude*rand(grid.size);
        end

        % create cell arrays with the same size as the elevation grid to
        % store cell connectivity (grid.flowsTo and grid.flowsFrom)
        grid.flowsTo = cell(grid.size);
        grid.flowsFrom = cell(grid.size);
        grid.flowsToCount = zeros(grid.size);
        
        % create grids to flag other attributes
        grid.oceanFlag = grid.z<=oceanLevel.z(1); % flag to identify whether a cell is in the ocean, using initial ocean level.
        grid.channelFlag = false(grid.size); % flag to identify whether a cell is part of a channel
    end
    