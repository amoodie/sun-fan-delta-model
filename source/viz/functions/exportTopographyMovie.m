function exportTopographyMovie(options)
% exportTopographyMovie.m: function that takes input data from the 
% Sun et al. (2002) fan-delta model and exports a movie of the topography.

fprintf('Starting movie export for %s ...\n',options.movieFilename)

% Load model parameters
load(options.modelParameterFile,'parameters');
runName = options.runName;
oceanLevel = parameters.oceanLevel; % meters
tSaveInterval_yr = parameters.tSaveInterval_yr;
tMax_yr = parameters.tMax_yr;
timesWithData_yr = 0:tSaveInterval_yr:tMax_yr;

% Load timeseries of DEMs and store in a structure array,
% dem.timeseries
for i=1:numel(timesWithData_yr)
    t = timesWithData_yr(i);
    filename = [options.basedir,'\',runName,'_time_',num2str(t,'%4.2f'),'_yr.mat'];
    load(filename,'grid')
    
    if i==1
        % save grids of x- and y- coordinates, which are constant across
        % the timeseries
        dem.x = grid.x;
        dem.y = grid.y;
    end
    dem.timeseries(i).z = grid.z;
    dem.timeseries(i).t = timesWithData_yr(i);
    dem.timeseries(i).hs = hillshade(dem.x,dem.y,dem.timeseries(i).z);
    
    % look up ocean level at this time
    indOceanLevel = find(dem.timeseries(i).t >=oceanLevel.timeStart_yr,1,'last');
    dem.timeseries(i).oceanLevel.z = oceanLevel.z(indOceanLevel);
end
clear grid

% Open video object
% See VideoWriter help for additional properties
v = VideoWriter([options.movieFilename,'.mp4'],'MPEG-4'); % Create video object v with the given save name
v.FrameRate = options.movieFrameRate; % Set the frames per second of the video
v.Quality = 100;            % Set the video quality (%)
open(v)                     % Open video object to write

% Create figure and axes
h = figure('Visible','on');    
set_plot_dimensions(options.movieFrameWidth,options.movieFrameHeight);
set(h,'color','k')

if options.setMinElevZero
    % Determine minimum elevation across timesteps
    minElevAll = Inf;
    for i=1:numel(dem.timeseries)
        elevs = dem.timeseries(i).z(:);
        minElev = min(elevs);
        if minElev<minElevAll
            minElevAll=minElev;
        end
    end
        
    % subtract minimum elevation (across all DEMs) from all DEMs to zero out
    for i=1:numel(dem.timeseries)
        dem.timeseries(i).oceanLevel.z = dem.timeseries(i).oceanLevel.z-minElevAll;
    end
    % Apply same offset to ocean level
    dem.timeseries(i).oceanLevel.z = dem.timeseries(i).oceanLevel.z - minElevAll;
end

% identify range of elevations across timesteps so can fix z-axis
% limits.
minElevAll = Inf;
maxElevAll = -Inf;
for i=1:numel(dem.timeseries)
    elevs = dem.timeseries(i).z(:);
    minElev = min(elevs);
    if minElev<minElevAll
        minElevAll=minElev;
    end
    maxElev = max(elevs);
    if maxElev>maxElevAll
        maxElevAll = maxElev;
    end
end
    
% Loop over topography timesteps and write video frames:
for i = 1:numel(dem.timeseries)
    % Subplot 1: 
    % - surf plot of *initial surface* in black 
    % - surf plot of shaded relief
    % - Plot transect location
    subplot(121)
    cla
    
    % Plot shaded relief
    s1=surf(dem.x/1000,dem.y/1000,dem.timeseries(i).z,dem.timeseries(i).hs);
    cax = [min(dem.timeseries(i).hs(:)) max(dem.timeseries(i).hs(:))];
    shading flat
    colormap gray
    
    if ~isnan(dem.timeseries(i).oceanLevel.z)
        % Plot ocean so that it extends way from surface
        oceanSurf = dem.timeseries(i).oceanLevel.z*ones(size(dem.x)); 
        oceanSurf(dem.timeseries(i).z> dem.timeseries(i).oceanLevel.z) = NaN;
        hold on, s2=surf(dem.x/1000,dem.y/1000,oceanSurf,'facecolor',[0 0 1],'facealpha',0.75,'edgecolor','none');
        % restore caxis for hillshade surface
        caxis(cax)
    end
    
    % set perspective view
    view([-60 10])
    % Axis labels
    xlabel('Distance (km)')
    ylabel('Distance (km)','fontsize',10)
    zlabel('Elevation (m)')
        
    % Profile
    if isequal(options.profileLimits.x(1),options.profileLimits.x(2))
        profile.x = repmat(options.profileLimits.x(1),1000,1);
    else 
        profile.x = linspace(options.profileLimits.x(1),options.profileLimits.x(2),1000)'; % i.e., 1000 points in profile
    end
    
    if isequal(options.profileLimits.y(1),options.profileLimits.y(2))
        profile.y = repmat(options.profileLimits.y(1),1000,1);
    else
        profile.y = linspace(options.profileLimits.y(1),options.profileLimits.y(2),1000)';
    end
    
    profile.z = interp2(dem.x,dem.y,dem.timeseries(i).z,profile.x,profile.y,'nearest'); 
    profile.d = [0;cumsum(sqrt(sum(diff([profile.x(:),profile.y(:)],1,1).^2,2)))]; % profile distance
    
    % Create a version of elevation profile with NaN where profile falls below ocean level
    ind =  profile.z < dem.timeseries(i).oceanLevel.z;
    profile.zClipOcean = profile.z;
    profile.zClipOcean(ind) = NaN;
    %hold on, plot3(profile.x/1000,profile.y/1000,profile.zClipOcean,'r-','linewidth',2)
    hold on, plot3(profile.x/1000,profile.y/1000,profile.z,'r-','linewidth',2)
    
    % create ghost handles for legend
    hold on, gh1(1) = plot(NaN,NaN,'linewidth',2,'color','r');
    hold on, gh1(2) = patch(NaN,NaN,'k');
    set(gh1(2),'facecolor',0.5*ones(1,3));
    hold on, gh1(3) = patch(NaN,NaN,'k');
    set(gh1(3),'facecolor','b','edgecolor','w','facealpha',0.75)
    
    if ~isnan(dem.timeseries(i).oceanLevel.z)
        legend(gh1(1:3),'\color{white}Profile','\color{white}Surface','\color{white}Standing water',...
            'Location','northeast','color',[0 0 0],'edgecolor',[0 0 0])
    else
        legend(gh1(1:2),'\color{white}Profile location','\color{white}Surface',...
            'Location','northeast','color',[0 0 0],'edgecolor',[0 0 0])
    end
        
    % set z-axis limits for consistency across timesteps
    set(gca,'zlim',[minElevAll maxElevAll])
    
    % set color to black, axis text color to white
    set(gca,'color','k')
    set(gca,'zcolor','w')
    set(gca,'ycolor','w')
    set(gca,'xcolor','w')
    set(gca,'fontsize',12)
    
    % Turn on plot box and grid to enhance three-dimensionality
    %grid on % incompatible with variable 'grid'
    box on
    
    % Reverse y-tick labels for consistency with second subplot
    if i==1
        yTickLabels = get(gca,'yticklabels');
        yTickLabels = flipud(yTickLabels);
    end
    set(gca,'yticklabels',yTickLabels)
    
    subplot(122)
    % Subplot 2:
    % - line plot of *initial surface* along profile (black)
    cla
    clear gh1 % clears a figure handle
    if i==1    
        profile.z_initial = interp2(dem.x,dem.y,dem.timeseries(1).z,profile.x,profile.y,'bilinear');
   %     profile.z_initial(profile.z_initial<dem.timeseries(i).oceanLevel.z) = NaN; 
    end
    
   plot(profile.d/1000,profile.z_initial,'-','linewidth',1,'color','w'); 
   %hold on, plot(profile.d/1000,profile.zClipOcean,'r-','linewidth',1)
   hold on, plot(profile.d/1000,profile.z,'r-','linewidth',1)
   
   
   bufferY = 0.1; % buffer for y-axis limits
    if ~isnan(dem.timeseries(i).oceanLevel.z)
        % - line plot of ocean level. Clip to intersection point with surface
        % profile.
        [~,~,io,~] = intersections(profile.d,dem.timeseries(i).oceanLevel.z*ones(size(profile.d)),profile.d,profile.z);
        if isempty(io) % i.e., no intersections  
            if all(profile.z<dem.timeseries(i).oceanLevel.z) % all topography is submerged
                oceanPatch.d = profile.d([1:end,end:-1:1]);
                oceanPatch.z = [dem.timeseries(i).oceanLevel.z*ones(size(profile.d));profile.z(end:-1:1)];
            elseif all(profile.z>dem.timeseries(i).oceanLevel.z)
                % Set ocean patch coordinates to empty as all topography is
                % subaerial
                 oceanPatch.d = [];
                 oceanPatch.z = [];
            end 
        else
            io = round(io);
            oceanPatch.d = [profile.d(io:end);profile.d(end:-1:io)];
            oceanPatch.z = [dem.timeseries(i).oceanLevel.z*ones(size(profile.d(io:end)));profile.z(end:-1:io)];
        end
        
        p = patch(oceanPatch.d/1000,oceanPatch.z,'k');
        set(p,'facecolor','b','edgecolor','none','facealpha',0.75)
    end
    xlabel('Distance (km)','fontsize',10)
    ylabel('Elevation (m)')
   
    % set x-axis limits
    set(gca,'xlim',[0 range(dem.y(:))/1000]) % y is the alongstream direction, which is the x-axis
    
    % set y-axis limits for consistency across timesteps
    set(gca,'ylim',[minElevAll-bufferY maxElevAll+bufferY])
    
    % create ghost handles for legend
    hold on, gh(1) = plot(NaN,NaN,'linewidth',2,'color','r');
    hold on, gh(2) = plot(NaN,NaN,'linewidth',1,'color','w');
    hold on, gh(3) = patch(NaN,NaN,'k');
    set(gh(3),'facecolor','b','edgecolor','w','facealpha',0.75)
    
    if ~isnan(dem.timeseries(i).oceanLevel.z)
        legend(gh(1:3),'\color{white} Profile','\color{white} Original profile','\color{white} Standing water',...
            'Location','northeast','color',[0 0 0],'edgecolor',[0 0 0])
    else
        legend(gh(1:2),'\color{white}Profile','\color{white}Original profile',...
            'Location','northeast','color',[0 0 0],'edgecolor',[0 0 0])
    end
    title(['\it{t}\rm = ',num2str(dem.timeseries(i).t,'%1.2f'),' yr'],'color','w')
     
    % set color to black, axis text color to white
    set(gca,'color','k')
    set(gca,'ycolor','w')
    set(gca,'xcolor','w')
    set(gca,'fontsize',12)
    
    % Add figure as a new frame in the video
    F = getframe(h);
    try
        writeVideo(v,F)
    catch
        error('Error writing video frame')
    end
    
   set(gcf,'invertHardCopy','off')
    % Also export movie frame as image
    print(gcf,[options.frameBasename,num2str(i,'%03.0f')],'-dpng')
end

% 7. Close figures and video object
close(v)
close(h)
close all
fprintf('Finished processing movie\n')
end