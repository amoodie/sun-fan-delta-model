% wrapper_run1.m: wrapper script to set parameters for run 1 and execute
% the model run using sunFanModel.m.

% we can load a checkpoint file and continue the model run
% to load a checkpoint, use a string pointing to a .mat filename output
% from a previous model run. To start a new run, use `false`.
loadCheckpoint = false;  % false | filename

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Start of parameters to edit
runName = 'sloped_run3'; % base name for run and file output
clobber = true; % whether to overwrite output folder if exists

% add the model source folder to the path
codeDir = genpath(fullfile('.','..','source'));
addpath(codeDir); % add source

% dimensional parameter constant set
con = loadConstants('mars-quartz-water');

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
Qw_inlet = 1000; % water discharge, m^3/s (in Table 2, base case: Qw_inlet = 20)
Qs_inlet = 2; % sediment discharge, m^3/s (named Q_sf in original paper. In Table 2, base case: 0.04)
Qw_threshold = 0.0000000005; % water discharge fraction to cut off channels
Qw_mismatch_tolerance = 1e-3; % tolerance param for raising a water mass-conservation error
Qs_threshold = Qs_inlet * 1e-9; % threshold amount of sediment transport for enacting an avulsion at cell
branchLimit = 8;

grid.dx = 500; % grid spacing, m (named "a" in the paper)
grid.xExtent = 200*grid.dx; % side length of square domain, m (named L_b in the paper)
grid.yExtent = grid.xExtent / 2; % added separate variabel for side length if y-dimension of grid
% Parameters for initial topography (can add additional options here)
grid.DEMoptions.initialSurfaceGeometry.type = 'slope'; % 'slopeBreak' | 'flat' % a 'flat' condition is used in Sun et al. (2002)
grid.DEMoptions.initialSurfaceGeometry.minElev = 0; 
grid.DEMoptions.addNoise = true;
grid.DEMoptions.noiseAmplitude = 0.001; % meters
grid.DEMoptions.slope.slope = -0.00083; % slope below slope break

% time paramaeters
t = 0; % Initial time
tMax_yr = 100; % simulation time, years
tStep_yr = 1e-4; % time step, years. Not specified in paper. There is some upper bound for stable topography change using the default input water/sediment discharges and grid cell spacing.
tSaveInterval_yr = 0.5; % time interval for saving data to file, years
tElapsedSinceSave_yr = 0; % variable to record elapsed time since saving

% boundary conditions
inlet.row = 1; % set inlet point for water and sediment
inlet.col = 100;
boundaryCondition = 'closed';
oceanLevel.steps = 500;
oceanLevel.timeStart_yr = linspace(0,tMax_yr,oceanLevel.steps); % times that define start of intervals with a particular ocean level
z0 = -1*grid.DEMoptions.slope.slope * grid.yExtent + grid.DEMoptions.initialSurfaceGeometry.minElev;
t0 = floor(oceanLevel.steps / 5) * 2;
zend = grid.DEMoptions.initialSurfaceGeometry.minElev;
tend = floor(oceanLevel.steps / 5) * 3;
oceanLevel.z = ones(oceanLevel.steps, 1);
oceanLevel.z(1:t0) = z0;
oceanLevel.z(tend+1:end) = zend;
oceanLevel.z(t0+1:tend) = linspace(z0, zend, floor(oceanLevel.steps / 5) * 1); % timeseries elevation of ponded water, m (xi_theta in the paper). The length of this vector must equal the length of the previous parameter.

% Flag to show a debugging figure. This is computationally expensive, so
% only use to debug.
debugFigure = false;

% set a rng seed for reproducible timing
rng(1)

%%% End of parameters to edit 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define name for output directory
outputDir = fullfile('.','..','output',runName);

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

% Execute model run, pausing if an error is thrown
dbstop if error
sunFanModel(parameters)
