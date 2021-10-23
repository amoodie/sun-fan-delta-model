function sunFanModel(parameters)
% sunFanModel.m: Implements numerical model in Sun et al. (2002), Fluvial 
% fan deltas: Linking channel processes with large-scale morphodynamics, 
% Water Resources Research 38(8), doi:10.1029/2001WR000284. All model
% parameters are taken from Tables 1 and 2 in the paper.
% Created June 30, 2021
% Ajay B. Limaye, University of Virginia (ajay@virginia.edu) 
% Andrew Moodie, UT Austin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process model parameters 

% Check that the parameters use the expected data types
checkParameters(parameters);

% extract parameters from structure array
runName = parameters.runName; 
clobber = parameters.clobber;
alpha_so = parameters.alpha_so;
alpha_sa = parameters.alpha_sa;
alpha_r = parameters.alpha_r;
alpha_b = parameters.alpha_b;
tauStar_c = parameters.tauStar_c;
n = parameters.n;
p = parameters.p;
R = parameters.R; 
g = parameters.g;
gamma = parameters.gamma;
lambda = parameters.lambda;
beta = parameters.beta;
Qw_inlet = parameters.Qw_inlet;
Qs_inlet = parameters.Qs_inlet;
Qw_threshold = parameters.Qw_threshold;
Qs_threshold = parameters.Qs_threshold;
Qw_mismatch_tolerance = parameters.Qw_mismatch_tolerance; 
D = parameters.D; 
oceanLevel = parameters.oceanLevel;
grid = parameters.grid;
t = parameters.t;
tStep_sec = parameters.tStep_sec;
tMax_yr = parameters.tMax_yr; 
tMax_sec = parameters.tMax_sec;
startingTime = parameters.startingTime;
tSaveInterval_yr = parameters.tSaveInterval_yr;
secondsPerYear = parameters.secondsPerYear;
tElapsedSinceSave_yr = parameters.tElapsedSinceSave_yr;
inlet = parameters.inlet;
boundaryCondition = parameters.boundaryCondition;
debugFigure = parameters.debugFigure;
loadCheckpoint = parameters.loadCheckpoint;
outputDir = parameters.outputDir;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set up an output folder for this run
if ~exist(outputDir,'dir')
       mkdir(outputDir)
else
    % error if clobber is set to false, so not to overwrite runs
    if ~clobber
        error('Output directory already exists and clobber is false')
    end
end

% save parameters to file
filename = fullfile(outputDir,[runName,'_parameters.mat']);
save(filename,'parameters')

% Create grid for elevation and cell connectedness, and grids of flags for channel cells, ocean cells
grid =  makeGrids(grid,oceanLevel);

% update slope, stored in 'grid'
grid = updateSlope(grid,boundaryCondition,t);

% flag grid cell with inlet as a channel cell
grid.channelFlag(inlet.row,inlet.col) = true;
grid.inletCell = sub2ind(size(grid.x),inlet.row,inlet.col);

% create grid to track elevation change in the last time step
grid.deltaz = zeros(grid.size); 

% handle loading the checkpoint, or creating the new output files
if ischar(loadCheckpoint)
    % load grid, oceanLevel, and startingTime from the checkpoint file
    data = load(loadCheckpoint);
    grid = data.grid;
    oceanLevel = data.oceanLevel;
    startingTime = data.t;
    clear data
else
    % save initial grid
    filename = fullfile(outputDir,[runName,'_time_0.00_yr.mat']);
    save(filename,'grid','t')
    fprintf('Saved file %s\n',filename)
end

% create a debugging figure if specified
if debugFigure
    debugFigureUpdateTime = tStep_sec * 100; % how often to update
    debugFig = figure('Position', [10 10 900 600]);
end

%%%%%%%% time loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iter = 0;
    for t = startingTime:tStep_sec:tMax_sec
        
        % compute the flowsToCount for each cell
        grid = countFlowsToInds(grid);
        
        %%% ARE THE ARRAYS THE SAME??
%         cnt = sum(grid.flowsToGraph, 1);
%         cnt = cnt(grid.channelFlag);
%         eql = grid.flowsToCount(~isnan(grid.flowsToCount)) == cnt;
%         if ~all(eql)
%             keyboard()
%         end

        % route flow to get discharge along each channel
        grid=routeFlow(grid,inlet,Qw_inlet,gamma,Qw_mismatch_tolerance);

        % look up slope along flow paths
        grid=slopeAlongFlow(grid); % updates field grid.S.alongFlow
        
        % update morphodynamic variables for each cell in the grid
        % (channel width, chanenl depth and sediment flux)
        grid = updateMorphodynamicVariables(grid,alpha_b,alpha_r,alpha_sa,alpha_so,R,g,D,tauStar_c,n,p,false);
        
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
        t_yr = t / secondsPerYear;
        indOceanLevel = find(t_yr>=oceanLevel.timeStart_yr,1,'last');
        grid.oceanLevel = oceanLevel.z(indOceanLevel);
        grid.oceanFlag = grid.z <= grid.oceanLevel;
        
        % check that any channels that are receiving flow below threshold
        % are disconnected from the network
        grid = unmarkAbandonedChannels(grid,Qw_threshold);

        % compute the flowsToCount for each cell (there may have been
        % changes to the network in above step)
        grid = countFlowsToInds(grid);
        
        % check for avulsion sites (criterion: eqn. 13). Change of flow
        % path from i-->j to i-->k initiated if criterion is met.
        avulsionCellInds = avulsionCheck(grid,beta,Qs_threshold);
        
        % if any avulsion sites were identified, enact avulsions that 
        % create new flow paths
        if ~isempty(avulsionCellInds)   
            grid = enactAvulsions(avulsionCellInds,grid,inlet);
        end
        
        % episodically save model output
        tElapsedSinceSave_yr = tElapsedSinceSave_yr + tStep_sec/(secondsPerYear);
        if tElapsedSinceSave_yr >= tSaveInterval_yr || t == tMax_yr
           filename = fullfile(outputDir,[runName,'_time_',num2str(t/(secondsPerYear),'%4.2f'),'_yr.mat']);
           save(filename,'grid','t')
           tElapsedSinceSave_yr = 0; % reset elapsed time since save 
           fprintf('Saved file %s\n',filename)
        end
        iter = iter + 1;
    end % end time loop
end % end function sunFanModel