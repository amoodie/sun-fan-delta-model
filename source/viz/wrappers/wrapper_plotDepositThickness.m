% wrapper_plotDepositThickness.m: wrapper script specifying plot options
% for plotDepositThickness.m.

clear
close all
dbstop if error

addpath ..\functions

% Specify 1 or 2 model cases to plot
options.modelCases.case1.runName = 'run1';
options.modelCases.case1.description = 'Fan-delta';
options.modelCases.case1.plotFlag = true;
options.modelCases.case1.dataDir = ['C:\Users\Ajay\Dropbox\work\active\development\Mars\Hypanis\sun-fan-delta-model\output\',options.modelCases.case1.runName,'\'];

options.modelCases.case2.runName = 'run1';
options.modelCases.case2.description = 'duplicate for test'; %'Alluvial fan';
options.modelCases.case2.plotFlag = true;
options.modelCases.case2.dataDir = ['C:\Users\Ajay\Dropbox\work\active\development\Mars\Hypanis\sun-fan-delta-model\output\',options.modelCases.case2.runName,'\'];

% set (x,y) limits for profile of sediment thickness
options.profileLimits.x = [5000 5000];
options.profileLimits.y = [10000 0]; 

% Specify output options
options.output.vizOutputDir = ['C:\Users\Ajay\Dropbox\work\active\development\Mars\Hypanis\sun-fan-delta-model\output\',options.modelCases.case1.runName,'\viz\'];
options.output.fileBasename = [options.output.vizOutputDir,options.modelCases.case1.runName,'_depositThickness'];
options.output.movie.frameRate = 2; % frames per second
options.output.movie.format = 'MPEG-4';
options.output.figureWidth = 20; % cm
options.output.figureHeight = 8; % cm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check that only 1 or 1&2 plots are defined
if ~options.modelCases.case1.plotFlag
    error('Define case 1 plot');
end

% Execute the plotting function using the specified options
if ~exist(options.output.vizOutputDir)
    mkdir(options.output.vizOutputDir);
end
plotDepositThickness(options);