function sunFanModel()
% sunFanModel.m: Implements numerical model in Sun et al. (2002), Fluvial 
% fan deltas: Linking channel processes with large-scale morphodynamics, 
% Water Resources Research 38(8), doi:10.1029/2001WR000284. All model
% parameters are taken from Tables 1 and 2 in the paper.
% Created June 30, 2021
% Ajay B. Limaye, University of Virginia (ajay@virginia.edu) 
% Andrew Moodie, UT Austin

clear,clc
dbstop if error

% add the model source folder to the path
addpath(genpath('source'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set model parameters 
runName = 'run1'; % base name for run and file output
clobber = false; % whether to overwrite output folder if exists

% Dimensionless parameters (from Table 1)
alpha_so = 11.25;
alpha_sa = 1;
alpha_r = 15;
alpha_b = 2.97;
tauStar_c = 0; %%% important parameter. Was set to 0 in Sun et al. (2002) because in
% that case deposition was always triggered by a standing water body.
% However, the paper assumes that tau* is always above critical (see eqn.
% 9). May need to modify this condition in some way to enable deposition -
% perhaps with a thin-flow cutoff for sediment transport (i.e., if depth is
% below a threshold value, set Qs to zero). 
n = 2.5;
p = 0;
R = 1.65; 
g = 3.7; % gravitational acceleration, m^2/s, for Mars. Gravity affects calculation of channel width (B) and sediment flux (Q_s). Do any of the other parameters depend on g?
gamma = 0.5;
lambda = 0.4;
beta = 1;

% Dimensioned parameters (from Table 2, base case)
Qw_inlet = 20; % water discharge, m^3/s
Qs_inlet = 0.04; % sediment discharge, m^3/s (named Q_sf in original paper)
Qw_threshold = 0.05; % water discharge fraction to cut off channels
Qw_mismatch_tolerance = 1e-3; % tolerance param for raising a water mass-conservation error
D = 0.3e-3; % grain diameter, m
oceanLevel = 0.01; %elevation of ponded water, m (named xi_theta in the paper)

grid.dx = 100; % grid spacing, m (named "a" in the paper)
grid.xExtent = 100*grid.dx; % side length of square domain, m (named L_b in the paper)
grid.yExtent = grid.xExtent; % added separate variabel for side length if y-dimension of grid

%%% Additional parameters. In the paper some of these are implicit or their values 
% are note given.
% time variables
t = 0; % Initialize time
tMax_yr = 10; % simulation time, years
tStep_yr = 1e-4; % time step, years. Not specified in paper. There is some upper bound for stable topography change using the default input water/sediment discharges and grid cell spacing.
tSaveInterval_yr = 0.1; % time interval for saving data to file, years
tElapsedSinceSave_yr = 0; % variable to record elapsed time since saving

% Parameters for initial topography
% start simple, will add additional options here
grid.DEMoptions.initialSurfaceGeometry.type = 'flat'; % 'slopeBreak' | 'flat' % a 'flat' condition is used in Sun et al. (2002)
grid.DEMoptions.initialSurfaceGeometry.minElev = 0; 
grid.DEMoptions.addNoise = true;
grid.DEMoptions.noiseAmplitude = 0.001; % meters

grid.DEMoptions.slopeBreak.fracDistToSlopeBreak = 0.5;
grid.DEMoptions.slopeBreak.carveChannel = true; 
grid.DEMoptions.slopeBreak.channelDepth = 2; % meters
grid.DEMoptions.slopeBreak.channelDepthTaperLength = 0.10; % fraction of array row dimension to use for tapering channel depth
grid.DEMoptions.slopeBreak.fracDistToSlopeBreak = 0.5; % Fractional distance along grid y-dimension at which slope break occurs 
grid.DEMoptions.slopeBreak.slope1 = -0.004; % slope above slope break
grid.DEMoptions.slopeBreak.slope2 = -0.001; % slope below slope break

% adjust tMax so that the value is reached in an integer number of time
% steps
tMax_yr = tStep_yr*ceil(tMax_yr/tStep_yr);

% Create grid for elevation and cell connectedness, and grids of flags for channel cells, ocean cells
grid =  makeGrids(grid,oceanLevel);

% set inlet point for water and sediment
% these values set the inlet to the low point in the first row of the DEM (i.e., the start of the channel)
inlet.row = 1; 
% [~,inlet.col] = min(grid.z(inlet.row,:));
inlet.col = 50;

boundaryCondition = 'closed';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grid = updateSlope(grid,boundaryCondition,t); % updates slope, stored in 'grid'

% flag grid cell with inlet as a channel cell
grid.channelFlag(inlet.row,inlet.col) = true;

grid.deltaz = zeros(grid.size); % change in elev per time

tStep_sec = tStep_yr*pi*1e7;
tMax_sec = tMax_yr*pi*1e7;

% set up an output folder for this run
outputFolderPath = fullfile('output',runName);
if ~exist(outputFolderPath,'dir')
       mkdir(outputFolderPath)
else
    % error if clobber is set to false, so not to overwrite runs
    if ~clobber
        error('Output folder already exists and clobber is false')
    end
end

% save initial grid
filename = fullfile(outputFolderPath,[runName,'_time_0.mat']);
save(filename,'grid','t')
fprintf('Saved file %s\n',filename)

% save model parameters (all the variables so far) to file
filename = fullfile(outputFolderPath,[runName,'_parameters.mat']);
save(filename)

% show a debugging figure?
% this is computationally expensive, so only use to debug
debugFigure = true;
debugFigureUpdateTime = tStep_sec * 10; % how often to update
if debugFigure
    debugFig = figure('Position', [10 10 900 600]);
end

%%%%%%%% time loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for t = tStep_sec:tStep_sec:tMax_sec
        
        % route flow to get discharge along each channel
        grid=routeFlow(grid,inlet,Qw_inlet,gamma,Qw_mismatch_tolerance);
        
        % look up slope along flow paths
        grid=slopeAlongFlow(grid); % updates field grid.S.alongFlow
        
        % calculate channel width (eqn. 9a), depth (eqn. 9b), sediment discharge
        % (eqn. 9c). Each of these is written in dimensionless form in the
        % corresponding equations, but most of their results plots (Fig. 2 to 7)
        % use dimensions. So here I calculate width, depth, and sediment discharge
        % in their dimensioned forms. Each quantity is calculated
        % simulataneously for all cells on the grid.
        
        % slope is defined as positive downhill.
        grid.B = (alpha_b^(-(3+2*p)/2))*(alpha_r^-1)*(grid.S.alongFlow.^(1+p))*((grid.Qw/(sqrt(g*D)*D^2)))*D; % width
        grid.H = (alpha_b./grid.S.alongFlow)*D; % depth
        % Note that given these dependencies on slope, both depth and width
        % can come out with non-physical values (i.e., negative or zero).
        % Those cases would also result in Qs_out = 0 based on the code
        % immediately below. The only other effect is that depth is used in
        % the avulsion check;  I think (but haven't double-check) cells with negative depth sd
        % cannot meet the avulsion criteria.
        
        grid.Qs_out = alpha_so*alpha_sa*(alpha_b^(-(3+2*p)/2))*(alpha_r^-1)*((alpha_b/R - tauStar_c)^n)*(grid.S.alongFlow.^(1+p)).*(grid.Qw/(sqrt(g*D)*D^2))*(sqrt(R*g*D)*D^2); % sediment discharge

        % set any points with Qs_out < 0 to Qs_out = 0 (can happen for
        % negative slope values (i.e., adverse slope)). In those cases,
        % no sediment should leave the cell.
        grid.Qs_out(grid.Qs_out<0) = 0;
        
        % update the debugging figure
        if debugFigure
            if mod(t, debugFigureUpdateTime) == 0
                debugFig = updateDebugFigure(debugFig, grid);
            end
        end
        
        % Enforce no sediment flux out for any cells that do not flow to
        % other cells
        grid.Qs_out(cellfun(@isempty,grid.flowsTo)) = 0;
        
        % Enforce no sediment flux out the for any cells below sea
        % level. (This is likely redundant, as no ocean cells should flow
        % to any other cells). 
        grid.Qs_out(grid.oceanFlag) = 0;
        
        if any(isnan(grid.Qs_out(grid.channelFlag)))
            error('Unexpected NaN value in sediment flux calculation');
        end
        
        % update bed elevation along flow paths defined by channel segments (eqn. 12)
        grid.Qs_in = zeros(grid.size);
        for i=1:grid.size(1)
            for j=1:grid.size(2)
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
                        for k=1:numel(contributingCells)
                            % Sediment routing at network junctions is not explictly
                            % discussed in the paper. For simplicitly, I assume
                            % that the sediment discharge is partitioned in the
                            % same way that water discharge is partitioned.

                            % look up Qs_out at the contributing cell; use the water discharge partition fraction at
                            % at the tributary to set the partition fraction
                            % for sediment
                            Qs_out_fromContributingCell = grid.Qs_out(contributingCells(k));
                            % Qs_out_fromContributingCell is sometimes empty:
                            % if discharge did not reach this location due to routeFlow rules,
                            % depsite these cells being connected by a channel in the tracking cell arrays.
                            % We thus include a safety here, to just proceed with the Qs calculation,
                            % but with 0 sediment contribution from this cell
                            if Qs_out_fromContributingCell > 0
                                Qw_fraction_receivedFromContributingCell = grid.flowsToFrac_Qw_distributed{contributingCells(k)};  
                                [receivingCellRow,receivingCellCol] = ind2sub(grid.size,grid.flowsTo{contributingCells(k)});
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
        end
        
        % set Qs_in, Qs_out to zero for non-channel cells (values are
        % currently NaN)
        grid.Qs_in(~grid.channelFlag)=0;
        grid.Qs_out(~grid.channelFlag)=0;
        
        % conservation of mass (equation 12)
        grid.deltaz = ((grid.Qs_in-grid.Qs_out)/((1-lambda)*grid.dx^2))*tStep_sec; % equivalent to  -((Qs_out-Qs_in)/((1-lambda)*grid.dx^2))*tStep; 
        
        % update topography for all grid cells
        grid.z = grid.z + grid.deltaz;
        
        grid = updateSlope(grid,boundaryCondition,t); % update slope arrays following topography update

        % update grid.oceanFlag following topography update
        grid.oceanFlag = grid.z <= oceanLevel;
        
        % check that any channels that are receiving flow below threshold
        % are disconnected from the network
        grid = unmarkAbandonedChannels(grid,Qw_threshold);
        
        % check for avulsion sites (criterion: eqn. 13). Change of flow
        % path from i-->j to i-->k initiated if criterion is met.
        
        % for each channel cell, look up the local bed elevation,
        % channel depth, and downstream slope. Compare to elevation for
        % neighboring cells and distance to those cells. Flag for avulsion
        % if criterion is met. I imposed an additional constraint that a
        % cell can only flow into 2 additional cells, so if the number of cells indicated grid.flowsTo is
        % already 2 or more, no new avulsions can be made.
        newAvulsions.rSource = []; % initialize structure array to store row and column coordinates for new avulsion cells
        newAvulsions.cSource = [];
        newAvulsions.indSource = [];
        newAvulsions.rNew = [];
        newAvulsions.cNew = [];
        newAvulsions.indNew = [];
        
        for k=1:numel(grid.channelFlag)
            if grid.channelFlag(k) && numel(grid.flowsTo{k})<2 % i.e., if it's a channel cell and flows to no more than 1 cell, then eligible for a new avulsion
                z_i = grid.z(k);
                H_ij = grid.H(k); % listed as "H_ij" in the paper, but I don't understand that notation as it seems to imply depth from i to j; easier to conceive as depth at i

                if isinf(H_ij)|| isnan(H_ij) || H_ij <= 0 
                    % because depth is calculated using slope, and slope can be negative, the value of 
                    % depth can be non-physical (Inf / NaN / negative). In that
                    % case, it cannot be used in the avulsion criterion in
                    % equation 13; I think a reasonable workaround is to
                    % set H_ij to zero for this case so that flow depth
                    % does not impact avulsion susceptibilty for this case.
                    H_ij = 0;
                end
                    
                S_ij = grid.S.alongFlow(k);
                % search nearest neighbor cells for potential avulsion
                % destinations. Skip any current channel cells or ocean
                % cells.
                [currentRow,currentCol] = ind2sub(grid.size,k);

                % search the 8 nearest neighbor cells
                rowSearch = currentRow + [-1 -1 -1 0 1 1 1 0]; % clockwise from NW
                colSearch = currentCol + [-1 0 1 1 1 0 -1 -1];
                L_ik = grid.dx*[sqrt(2) 1 sqrt(2) 1 sqrt(2) 1 sqrt(2) 1]; % distance from starting cell to search cell
                avulsionSusceptibilityIndex = nan(1,8); % NaN rather than 0 so that "0" doesn't accidentally become a maximum
                for l=1:numel(rowSearch)
                    if rowSearch(l) < 1 || rowSearch(l) > grid.size(1) || colSearch(l)<1 || colSearch(l) > grid.size(2) || grid.channelFlag(rowSearch(l),colSearch(l))
                        continue % if search cell is already a channel cell or is off the grid then cannot avulse to there. Therefore, leave avulsion susceptibility index as NaN for this cell and continue to next loop iteration.
                    else
                        z_k = grid.z(rowSearch(l),colSearch(l));
                        avulsionSusceptibilityIndex(l) = ((z_i-beta*H_ij) - z_k)/L_ik(l); % equation 13 LHS
                    end
                end

                if any(avulsionSusceptibilityIndex > S_ij) % equation 13
                    % select the neighboring cell with the greatest avulsion
                    % susceptibility 
                    [~,neighborAvulsionSelect] = max(avulsionSusceptibilityIndex);
                    newAvulsions.rNew = [newAvulsions.rNew; rowSearch(neighborAvulsionSelect)];
                    newAvulsions.cNew = [newAvulsions.cNew; colSearch(neighborAvulsionSelect)];
                    indNew = sub2ind(grid.size,rowSearch(neighborAvulsionSelect),colSearch(neighborAvulsionSelect));
                    newAvulsions.indNew = [newAvulsions.indNew; indNew];
                    newAvulsions.rSource = [newAvulsions.rSource; currentRow];
                    newAvulsions.cSource = [newAvulsions.cSource; currentCol];
                    newAvulsions.indSource = [newAvulsions.indSource; k];
                end
            end
        end
        % generate a new channel path originating from each avulsion site 
        % (paragraph 20 and eqn. 15). The path construction stops when any of the
        % following conditions are met: (1) path meets an existing channel, 
        % (2) path reaches standing water (a "wet" cell), or (3) path runs into a
        % topographic sink (local minimum). 
        
        if ~isempty(newAvulsions.rNew)
            
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
        end % end path construction for avulsions
        
        % episodically save model output
        tElapsedSinceSave_yr = tElapsedSinceSave_yr + tStep_sec/(pi*1e7);
        if tElapsedSinceSave_yr >= tSaveInterval_yr || t == tMax_yr
           filename = fullfile(outputFolderPath,[runName,'_time_',num2str(t/(pi*1e7),'%4.2f'),'_yr.mat']);
           save(filename,'grid','t','oceanLevel')
           tElapsedSinceSave_yr = 0; % reset elapsed time since save 
           fprintf('Saved file %s\n',filename)
        end
    end % end time loop
end % end function sunFanModel