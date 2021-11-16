function grid = updateSlope(grid,boundaryCondition,t) % nested function
        % updateSlope: calculates slope values. Convention: slopes are
        % positive toward lower elevations. 

        grid.zPadded = grid.z;
        switch boundaryCondition
            case 'closed'
                grid.zPadded = padarray(grid.z,[1 1],NaN,'both'); % pads array with 1 NaN value
            case 'periodicXY'
                grid.zPadded = [grid.zPadded(1,:);grid.zPadded;grid.zPadded(end,:)];
                grid.zPadded = [grid.zPadded(:,1),grid.zPadded,grid.zPadded(:,end)];
        end

        %calculate slopes
        grid.S.d8(:) = NaN;
        grid.S.d8(1, :, :) = (grid.z - grid.zPadded(1:end-2, 1:end-2)) / (grid.dx * sqrt(2)); % NE
        grid.S.d8(2, :, :) = (grid.z - grid.zPadded(1:end-2, 2:end-1)) / (grid.dx);
        grid.S.d8(3, :, :) = (grid.z - grid.zPadded(1:end-2, 3:end)) / (grid.dx * sqrt(2));
        grid.S.d8(4, :, :) = (grid.z - grid.zPadded(2:end-1, 3:end)) / grid.dx;
        grid.S.d8(5, :, :) = (grid.z - grid.zPadded(3:end, 3:end)) / (grid.dx * sqrt(2));
        grid.S.d8(6, :, :) = (grid.z - grid.zPadded(3:end, 2:end-1)) / grid.dx;
        grid.S.d8(7, :, :) = (grid.z - grid.zPadded(3:end, 1:end-2)) / (grid.dx * sqrt(2));
        grid.S.d8(8, :, :) = (grid.z - grid.zPadded(2:end-1, 1:end-2)) / (grid.dx);

    end
