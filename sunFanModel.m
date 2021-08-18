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

% we can load a checkpoint file and continue the model run
% to load a checkpoint, use a string pointing to a .mat filename output
% from a previous model run. To start a new run, use `false`.
loadCheckpoint = false; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set model parameters 
runName = 'run5'; % base name for run and file output
clobber = true; % whether to overwrite output folder if exists

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
startingTime = tStep_sec;
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

% handle loading the checkpoint, or creating the new output files
if ischar(loadCheckpoint)
    % load from file the checkpoint
    data = load('output/run5/run5_time_0.40_yr.mat');
    grid = data.grid;
    oceanLevel = data.oceanLevel;
    startingTime = data.t;
else
    % save initial grid^M
    filename = fullfile(outputFolderPath,[runName,'_time_0.00.mat']);
    save(filename,'grid','t')
    fprintf('Saved file %s\n',filename)

    % save model parameters (all the variables so far) to file
    filename = fullfile(outputFolderPath,[runName,'_parameters.mat']);
    save(filename)
end

% show a debugging figure?
% this is computationally expensive, so only use to debug
debugFigure = true;
debugFigureUpdateTime = tStep_sec * 10; % how often to update
if debugFigure
    debugFig = figure('Position', [10 10 900 600]);
end

%%%%%%%% time loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for t = startingTime:tStep_sec:tMax_sec
        
        % check that the network is valid (and repair/trim) what is not
        grid=validateNetwork(grid,inlet);

        % route flow to get discharge along each channel
        grid=routeFlow(grid,inlet,Qw_inlet,gamma,Qw_mismatch_tolerance);

        % look up slope along flow paths
        grid=slopeAlongFlow(grid); % updates field grid.S.alongFlow
        
        % update morphodynamic variables for each cell in the grid
        % (channel width, chanenl depth and sediment flux)
        grid = updateMorphodynamicVariables(grid,alpha_b,alpha_r,alpha_sa,alpha_so,R,g,D,tauStar_c,n,p);
        
        % update the debugging figure
        if debugFigure
            if mod(t, debugFigureUpdateTime) == 0
                debugFig = updateDebugFigure(debugFig, grid);
            end
        end
        
        % update bed elevation along flow paths
        grid = updateTopography(grid,inlet,lambda,tStep_sec,Qs_inlet); 
        
        % update slope arrays following topography update
        grid = updateSlope(grid,boundaryCondition,t); 

        % update grid.oceanFlag following topography update
        grid.oceanFlag = grid.z <= oceanLevel;
        
        % check that any channels that are receiving flow below threshold
        % are disconnected from the network
        grid = unmarkAbandonedChannels(grid,Qw_threshold);
        
        % check for avulsion sites (criterion: eqn. 13). Change of flow
        % path from i-->j to i-->k initiated if criterion is met.
        newAvulsions = avulsionCheck(grid,beta);
        
        % if any avulsion sites were identified, enact avulsions that 
        % create new flow paths
        if ~isempty(newAvulsions.rNew)   
            grid = enactAvulsions(newAvulsions,grid,inlet);
        end
        
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