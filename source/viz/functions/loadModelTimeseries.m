function dem = loadModelTimeseries(modelParameterFile,dataDir)
% loadModelTimeseries.m: Loads timeseries data for a model run specified by
% modelParameterFile.

load(modelParameterFile,'parameters');
tSaveInterval_yr = parameters.tSaveInterval_yr;
tMax_yr = parameters.tMax_yr;
timesWithData = 0:tSaveInterval_yr:tMax_yr;
runName = parameters.runName;
oceanLevel = parameters.oceanLevel;

% Load timeseries of DEMs and store in a structure array,
% dem.timeseries
for i=1:numel(timesWithData)
    filename = [dataDir,runName,'_time_',num2str(timesWithData(i),'%12.2f'),'_yr.mat'];
    load(filename,'grid')
    
    if i==1
        % save grids of x- and y- coordinates, which are constant across
        % the timeseries
        dem.x = grid.x;
        dem.y = grid.y;
        z_initial = grid.z;
    end
    dem.timeseries(i).z = grid.z;
    dem.timeseries(i).t = timesWithData(i);
    dem.timeseries(i).depositThickness = grid.z - z_initial;
    
    % look up ocean level at this time
    indOceanLevel = find( dem.timeseries(i).t >=oceanLevel.timeStart_yr,1,'last');
    dem.timeseries(i).oceanLevel.z = oceanLevel.z(indOceanLevel);
end

% get maximum deposit thickenss for the timesieres, which is used to set
% the caxis range for plotting
maxThickness = -Inf;
for i=1:numel(dem.timeseries)
    maxThicknessTemp = max(dem.timeseries(i).depositThickness(:));
    if maxThicknessTemp > maxThickness
        maxThickness = maxThicknessTemp;
    end
end
   
dem.maxSedThickness = maxThickness;

end