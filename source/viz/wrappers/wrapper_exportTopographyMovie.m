% wrapperExportMovie: Wrapper script for exporting movies of topography.

clear,close all,clc
dbstop if error

addpath ..\functions

% % Run 1
options.runName = 'run1';
options.setMinElevZero = true; % uses minimum elevation in the set of DEMs to define zero
options.movieFrameRate = 2; % frames per second for movie
options.movieFrameWidth = 24; % cm
options.movieFrameHeight = 8; % cm
options.importedDataFile = [options.runName,'data.mat'];
options.basedir = ['C:\Users\Ajay\Dropbox\work\active\development\Mars\Hypanis\sun-fan-delta-model\output\',options.runName];
options.modelParameterFile = [options.basedir,'\',options.runName,'_parameters.mat'];
options.vizOutputDir = [options.basedir,'\viz\'];
options.movieFilename = [options.vizOutputDir,options.runName,'_topoEvolution'];
options.frameBasename = [options.vizOutputDir,options.runName,'_topoEvolution_frame'];

% set (x,y) limits for profile of surface
options.profileLimits.x = [5000 5000];
options.profileLimits.y = [10000 0]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist(options.vizOutputDir)
    mkdir(options.vizOutputDir);
end
exportTopographyMovie(options);