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
        
        if t==0
            % Initialize stencils that will be used to select elevation/slope/discharge
            % values from padded versions of the grid
              [nrows,ncols] = size(grid.zPadded);
              grid.stencil.NW.r = 1:nrows-2;
              grid.stencil.NW.c = 1:ncols-2;
              grid.stencil.N.r = 1:nrows-2;
              grid.stencil.N.c = 2:ncols-1;
              grid.stencil.NE.r = 1:nrows-2;
              grid.stencil.NE.c = 3:ncols;
              grid.stencil.E.r = 2:nrows-1;
              grid.stencil.E.c = 3:ncols;
              grid.stencil.SE.r = 3:nrows;
              grid.stencil.SE.c =  3:ncols;
              grid.stencil.S.r = 3:nrows;
              grid.stencil.S.c = 2:ncols-1;
              grid.stencil.SW.r = 3:nrows;
              grid.stencil.SW.c = 1:ncols-2;
              grid.stencil.W.r = 2:nrows-1;
              grid.stencil.W.c = 1:ncols-2;
        end        
                
        %calculate slopes
        grid.S.NW = -(grid.zPadded(grid.stencil.NW.r,grid.stencil.NW.c)-grid.z)/(sqrt(2)*grid.dx);
        grid.S.N = -(grid.zPadded(grid.stencil.N.r,grid.stencil.N.c)-grid.z)/grid.dx;
        grid.S.NE = -(grid.zPadded(grid.stencil.NE.r,grid.stencil.NE.c)-grid.z)/(sqrt(2)*grid.dx);
        grid.S.E = -(grid.zPadded(grid.stencil.E.r,grid.stencil.E.c)-grid.z)/grid.dx;
        grid.S.SE = -(grid.zPadded(grid.stencil.SE.r,grid.stencil.SE.c)-grid.z)/(sqrt(2)*grid.dx);
        grid.S.S = -(grid.zPadded(grid.stencil.S.r,grid.stencil.S.c)-grid.z)/grid.dx;
        grid.S.SW = -(grid.zPadded(grid.stencil.SW.r,grid.stencil.SW.c)-grid.z)/(sqrt(2)*grid.dx);
        grid.S.W = -(grid.zPadded(grid.stencil.W.r,grid.stencil.W.c)-grid.z)/(sqrt(2)*grid.dx);
        
        switch boundaryCondition
            case 'closed'
                grid.S.NW(isnan(grid.S.NW)) = -Inf;
                grid.S.N(isnan(grid.S.N)) = -Inf;
                grid.S.NE(isnan(grid.S.NE)) = -Inf;
                grid.S.E(isnan(grid.S.E)) = -Inf;
                grid.S.SE(isnan(grid.S.SE)) = -Inf;
                grid.S.S(isnan(grid.S.S)) = -Inf;
                grid.S.SW(isnan(grid.S.SW)) = -Inf;
                grid.S.W(isnan(grid.S.W)) = -Inf;
        end
    end
