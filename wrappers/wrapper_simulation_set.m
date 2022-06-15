% wrapper_simulation_set.m: wrapper script to run a set of simulations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start of parameters to edit

setName = 'hypanis_simple_set0'; % base name for run and file output
clobber = true; % whether to overwrite output folder if exists

% add the model source folder to the path
codeDir = genpath(fullfile('.','..','source'));
addpath(codeDir); % add source

% dimensional parameter constant set
con = loadConstants('mars-basalt-water');

% Dimensionless parameters specified in the paper
alpha_so = 11.25; % from Table 1
alpha_sa = 1; % from Table 1
alpha_r = 15; % from Table 1
alpha_b = 2.97; % from Table 1
tauStar_c = 0; % from Table 1. This parameter is important and was set to 0 in Sun et al. (2002) because in
% that case deposition was always triggered by a standing water body.
% However, the paper assumes that tau* is always above critical (see eqn.
% 9). May need to modify this condition in some way to enable deposition -
% perhaps with a thin-flow cutoff for sediment transport (i.e., if depth is
% below a threshold value, set Qs to zero). 
n = 2.5; % from Table 1
p = 0; % from Table 1
R = con.R; % from Table 1
g = con.g; % gravitational acceleration, m^2/s, modified for Mars. Gravity affects calculation of channel width (B) and sediment flux (Q_s). Do any of the other parameters depend on g?
gamma = 0.5; % from Table 1
lambda = 0.4; % from Table 1
beta = 1; % from Table 1

% some dimensional params
D = 0.3e-3; % grain diameter, m (in Table 2, base case: D = 0.3e-3)

%%% Additional parameters. Some are implicit or unspecified in the paper;
%%% others are added for this model implementation. 

% Flow routing
Qw_inlet = 7500; % water discharge, m^3/s (in Table 2, base case: Qw_inlet = 20)
Qs_inlet = 10; % sediment discharge, m^3/s (named Q_sf in original paper. In Table 2, base case: 0.04)
Qw_threshold = 0.05; % water discharge fraction to cut off channels
Qw_mismatch_tolerance = 1e-3; % tolerance param for raising a water mass-conservation error
Qs_threshold = Qs_inlet * 0.05; % threshold amount of sediment transport for enacting an avulsion at cell
branchLimit = 2;

grid.dx = 1000; % grid spacing, m (named "a" in the paper)
grid.xExtent = 200*grid.dx; % side length of square domain, m (named L_b in the paper)
grid.yExtent = grid.xExtent / 2; % added separate variabel for side length if y-dimension of grid
% Parameters for initial topography (can add additional options here)
grid.DEMoptions.initialSurfaceGeometry.type = 'slope'; % 'slopeBreak' | 'flat' % a 'flat' condition is used in Sun et al. (2002)
grid.DEMoptions.initialSurfaceGeometry.minElev = 0; 
grid.DEMoptions.addNoise = true;
grid.DEMoptions.noiseAmplitude = 0.1; % meters
grid.DEMoptions.slope.slope = -0.00083; % slope below slope break

% time paramaeters
t = 0; % Initial time
tMax_yr = 400; % simulation time, years
tStep_yr = 1e-3; % time step, years. Not specified in paper. There is some upper bound for stable topography change using the default input water/sediment discharges and grid cell spacing.
tSaveInterval_yr = 1; % time interval for saving data to file, years
tElapsedSinceSave_yr = 0; % variable to record elapsed time since saving

% boundary conditions
inlet.row = 1; % set inlet point for water and sediment
inlet.col = 100;
boundaryCondition = 'closed';

% timeseries elevation of ponded water, m (xi_theta in the paper).
%   we set this up as an array of directives for the rate of fall, which
%   are then transformed to a nt x n_runs array of oceanLevel
setRunOceanLevels = [NaN, 0]; % no water, constant, n_rates; (m/yr)
%   configure oceanLevel.timeStart_yr and oceanLevel.z to be the same
%   length, and jointly defining the water level curve. In the model,
%   water level is discretized by this curve, so choose a sufficient number
%   of steps that the curve is well defined and changing at least a few
%   times per year.
oceanLevel.steps = 1000;
oceanLevel.timeStart_yr = linspace(0,tMax_yr,oceanLevel.steps); % times that define start of intervals with a particular ocean level
oceanLevel.z = NaN(oceanLevel.steps, 1);  % preallocate as NaN (replace before starting run!)
oceanLevel.fallRate = NaN; % preallocate as NaN (replace before starting run!)
ntChunk = 4; % how to break up the fall and flats
tChunk = floor(oceanLevel.steps / ntChunk);  % duration of each chunk
tStartFall = 3;  % when to start fall, measured in tChunks
t0 = tStartFall*tChunk;  % index to start the fall
tFallTime = ntChunk - tStartFall; % how logn the fall, measured in tChunks
%   now, set up the array of oceanLevels for each run (preallocate as NaN)
oceanLevelArray = NaN(oceanLevel.steps, length(setRunOceanLevels));
for ol_i = 1:length(setRunOceanLevels)
    % preallocate the z array for this run
    zi = NaN(oceanLevel.steps, 1);
    % us the value to determine the starting condition for water
    %     if the value is nan, there should be no water at all
    if isnan(setRunOceanLevels(ol_i))
        zStartWater = NaN;
    else
        zStartWater = -1*grid.DEMoptions.slope.slope * (grid.yExtent-grid.dx) + grid.DEMoptions.initialSurfaceGeometry.minElev;
    end
    % make the whole array the start value
    zi(:) = zStartWater;
    % if there is a rate, then step down at this rate
    if ~isnan(setRunOceanLevels(ol_i))
        % determine the falling values for the whole run duration
        fallRate = setRunOceanLevels(ol_i); % m/yr
        if isfinite(fallRate)
            zEndIfInf = zStartWater - ((fallRate * (tMax_yr / oceanLevel.steps)) * tChunk * tFallTime);  % if rate continued until end of run
            % fill the values in the array
            zi(t0+1:end) = linspace(zStartWater, zEndIfInf, tChunk * tFallTime);
            % and replace everything below
            %    minimum with NaN
            zi(zi < grid.DEMoptions.initialSurfaceGeometry.minElev) = NaN;
        else
            % special case for Inf, this is the instant drop to half
            zEndInstant = zStartWater / 2;
            zi(t0+1:end) = zEndInstant;
        end
    end
    % place the temp array into the storage array
    oceanLevelArray(:, ol_i) = zi;
end

% did this work correctly?
% figure()
% plot(oceanLevelArray)

% Flag to show a debugging figure. This is computationally expensive, so
% only use to debug.
debugFigure = false;

% set a rng seed for reproducible timing
% rng(1)

%%% End of parameters to edit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set up the parallel loop and run the simulations

% set up the output folder for the set
setOutputDir = fullfile('.','..','output',setName);
% set up an output folder for this run
if ~exist(setOutputDir,'dir')
       mkdir(setOutputDir)
end

% make a placeholder for output folder of run
outputDir = fullfile(setOutputDir,'placeholder');
loadCheckpoint = false;

% adjust time parameters for use in model execution
tMax_yr = tStep_yr*ceil(tMax_yr/tStep_yr);% adjust tMax so that the value is reached in an integer number of time steps
secondsPerYear = 60*60*24*365; % seconds in one calendar year
tStep_sec = tStep_yr*secondsPerYear; % time step in seconds
startingTime = tStep_sec; % starting time (overwritten if loads existing model data from checkpoint)
tMax_sec = tMax_yr*secondsPerYear; % maximum time in seconds

%%% Pack all the variables together into a structure array to simplify I/O.
parametersCell = who;
parametersCell{end+1}='fieldNames'; % Needed for v2struct.m to use the variable names as structure fields.
parameters = v2struct(parametersCell);
% trim the parameters field
parameters = rmfield(parameters, {'oceanLevelArray', 'setOutputDir'});

parfor run_i = 1:length(setRunOceanLevels)
    % make a string of the run name
    runName = ['run', num2str(run_i-1)];
    disp(runName)
    
    % make a copy of the parameters
    parameters_i = parameters
    
    % replace values in parameters_i for this run
    parameters_i.runName = runName;
    parameters_i.oceanLevel.fallRate = setRunOceanLevels(run_i);
    parameters_i.oceanLevel.z = oceanLevelArray(:, run_i);
    parameters_i.outputDir = fullfile(setOutputDir, runName)

    % Execute model run, pausing if an error is thrown
    dbstop if error
    sunFanModel(parameters_i)
end

% compile the results
parfor run_i = 1:length(setRunOceanLevels)
    % make a string of the run name
    runName = ['run', num2str(run_i-1)];
    disp(runName)
    
    % make a copy of the parameters
    compileResultsToNetCDF(fullfile(setOutputDir, runName))
end

