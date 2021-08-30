function plotDepositThickness(options)
% plotDepositThickness.m: plots deposit thickness for 1 or 2 model cases:
% - contours of sediment thickness
% - raster of sediment thickness
% - outputs: plots of final sediment thickness as well as movies of sediment thickness over time

% Create video objects
% See VideoWriter help for additional properties
v_sedimentThicknessContours = VideoWriter([options.output.fileBasename,'_contourTimeseries'],options.output.movie.format); % Create video object v with the given save name
v_sedimentThicknessContours.FrameRate = options.output.movie.frameRate; % Set the frames per second of the video
v_sedimentThicknessContours.Quality = 100; % Set the video quality (%)
open(v_sedimentThicknessContours) % Open video object to write

v_sedimentThicknessProfile = VideoWriter([options.output.fileBasename,'_profileTimeseries'],options.output.movie.format); % Create video object v with the given save name
v_sedimentThicknessProfile.FrameRate = options.output.movie.frameRate; % Set the frames per second of the video
v_sedimentThicknessProfile.Quality = 100; % Set the video quality (%)
open(v_sedimentThicknessProfile) % Open video object to write

v_sedimentThicknessRaster = VideoWriter([options.output.fileBasename,'_rasterTimeseries'],options.output.movie.format); % Create video object v with the given save name
v_sedimentThicknessRaster.FrameRate = options.output.movie.frameRate; % Set the frames per second of the video
v_sedimentThicknessRaster.Quality = 100; % Set the video quality (%)
open(v_sedimentThicknessRaster) % Open video object to write

% Create figures
fh_sedimentThicknessContours = figure('Visible','on');    
set_plot_dimensions(options.output.figureWidth,options.output.figureHeight);

% Create figures
fh_sedimentThicknessProfile= figure('Visible','on');    
set_plot_dimensions(options.output.figureWidth,options.output.figureHeight);

% Create figures
fh_sedimentThicknessRaster = figure('Visible','on');    
set_plot_dimensions(options.output.figureWidth,options.output.figureHeight);

% If plotting 2 model cases, set flag to create a 1x2 subplot
if options.modelCases.case1.plotFlag && options.modelCases.case2.plotFlag
    use2subplots = true;
else
    use2subplots = false;
end

% Load timeseries data for first model case
modelParameterFile = [options.modelCases.case1.dataDir,options.modelCases.case1.runName,'_parameters.mat'];
case1_data = loadModelTimeseries(modelParameterFile,options.modelCases.case1.dataDir);

if options.modelCases.case2.plotFlag
    % Load timeseries data for the (optional) second model case
    modelParameterFile = [options.modelCases.case2.dataDir,options.modelCases.case2.runName,'_parameters.mat'];
    case2_data  = loadModelTimeseries(modelParameterFile,options.modelCases.case2.dataDir);
end

% Plot sediment thickness using contours
figure(fh_sedimentThicknessContours)
for i=1:numel(case1_data.timeseries) % assumes both case1 and case2 have same number of timeseries entries
    timeseriesIndex = i;
    if use2subplots
        subplot(121)
        plotThicknessMapview(case1_data,timeseriesIndex,options.modelCases.case1.description,'contours');
        
        subplot(122)
        plotThicknessMapview(case2_data,timeseriesIndex,options.modelCases.case2.description,'contours');
    else
        plotThicknessMapview(case1_data,timeseriesIndex,options.modelCases.case1.description,'contours');
    end
    F = getframe(fh_sedimentThicknessContours);
    try
        writeVideo(v_sedimentThicknessContours,F)
    catch
        error('Error writing video frame')
    end
end
% export final image to separate file
print(gcf,[options.output.fileBasename,'_contoursFinal'],'-dpng')
close(fh_sedimentThicknessContours) % close figure
close(v_sedimentThicknessContours) % close video object

% Plot sediment thickness using rasters
figure(fh_sedimentThicknessRaster)
for i=1:numel(case1_data.timeseries) % assumes both case1 and case2 have same number of timeseries entries
    timeseriesIndex = i;
    if use2subplots
        subplot(121)
        plotThicknessMapview(case1_data,timeseriesIndex,options.modelCases.case1.description,'raster');
        
        subplot(122)
        plotThicknessMapview(case2_data,timeseriesIndex,options.modelCases.case2.description,'raster');
    else
        plotThicknessMapview(case1_data,timeseriesIndex,options.modelCases.case1.description,'raster');
    end
    F = getframe(fh_sedimentThicknessRaster);
    try
        writeVideo(v_sedimentThicknessRaster,F)
    catch
        error('Error writing video frame')
    end
end
% export final image to separate file
print(gcf,[options.output.fileBasename,'_rasterFinal'],'-dpng')
close(fh_sedimentThicknessRaster) % close figure
close(v_sedimentThicknessRaster) % close video object

% Plot sediment thickness profiles
figure(fh_sedimentThicknessProfile)
for i=1:numel(case1_data.timeseries) % assumes both case1 and case2 have same number of timeseries entries
    timeseriesIndex = i;
    if use2subplots
        subplot(121)
        plotThicknessProfile(case1_data,timeseriesIndex,options.profileLimits,options.modelCases.case1.description,'k'); 
        subplot(122)
        plotThicknessProfile(case2_data,timeseriesIndex,options.profileLimits,options.modelCases.case2.description,'k'); 
    else
        plotThicknessProfile(case1_data,timeseriesIndex,options.profileLimits,options.modelCases.case1.description,'k'); 
    end
    F = getframe(fh_sedimentThicknessProfile);
    try
        writeVideo(v_sedimentThicknessProfile,F)
    catch
        error('Error writing video frame')
    end
end

% Export image
print(gcf,[options.output.fileBasename,'_profileComparisonFinal'],'-dpng')

% Close figure and video 
close(fh_sedimentThicknessProfile)
close(v_sedimentThicknessProfile)

    function plotThicknessMapview(dem_data,timeseriesIndex,plotTitle,plotType)
        switch plotType
            case 'contours'
                contour(dem_data.x/1000,dem_data.y/1000,dem_data.timeseries(timeseriesIndex).depositThickness)
            case 'raster'
                imagesc(dem_data.x(:)/1000,dem_data.y(:)/1000,dem_data.timeseries(timeseriesIndex).depositThickness)
        end
        xlabel('Distance (km)')
        ylabel('Distance (km)')
        title(plotTitle)
        axis image
        set(gca,'xdir','normal')
        set(gca,'ydir','normal')
        colormap hot
        cbar = colorbar;
        ylabel(cbar,'Deposit thickness (m)')
        set(gca,'fontsize',10)
        caxis([0 dem_data.maxSedThickness])
        
        % could add: caxis, xlim, switching of x-direction coordinates (update
        % dem.x to be distance from slope break, for example)
    end

    function plotThicknessProfile(dem_data,timeseriesIndex,profileLimits,plotTitle,lineColor)
        xi = linspace(profileLimits.x(1),profileLimits.x(2),1000); % i.e., 1000 points in profile
        yi = linspace(profileLimits.y(1),profileLimits.y(2),1000);
        thickness = interp2(dem_data.x,dem_data.y,dem_data.timeseries(timeseriesIndex).depositThickness,xi,yi,'nearest'); 
        d = [0;cumsum(sqrt(sum(diff([xi(:),yi(:)],1,1).^2,2)))]; % profile distance
        plot(d/1000,thickness,'-','color',lineColor,'linewidth',1)
        xlabel('Distance (km)')
        ylabel('Deposit thickness (m)')
        set(gca,'fontsize',10)
        set(gca,'ylim',[0 dem_data.maxSedThickness])
        title(plotTitle)
    end
end