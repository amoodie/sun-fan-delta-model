function compileResultsToNetCDF(resultsFolder)
%% compile results into a single netcdf file timeseries.
% input should be the folder containing all of the results files
% creates a file called 'elevation_timeseries.nc' in the resultsFolder.

    %% define the file
    
    % collect all files in the folder
    outputFileList = dir(fullfile(resultsFolder,'*_time_*.mat'));

    % sort the list of names by the output time
    [~, reindex] = sort( str2double( regexp( {outputFileList.name}, '\d++', 'match', 'once' )));
    outputFileList = outputFileList(reindex);
    
    % open the first file to get the info for the netCDF file
    frst = load(fullfile(outputFileList(1).folder, outputFileList(1).name));
    xSize = numel(frst.grid.xVec);
    ySize = numel(frst.grid.yVec);
    
    % use the file index to determine the times that we will save
    timesToSave = str2double( regexp( {outputFileList.name}, '\d+\.\d', 'match', 'once' ));

    % set up the netcdf file
    netCDFpath = fullfile(resultsFolder, 'elevation_timeseries.nc');
    ncid = netcdf.create(netCDFpath,'NETCDF4');

    dim_length = netcdf.defDim(ncid,'length',ySize);
    dim_width = netcdf.defDim(ncid,'width',xSize);
    dim_time = netcdf.defDim(ncid,'time',length(outputFileList));
         
    colMajorDims = [dim_time, dim_length, dim_width];
    rowMajorDims = fliplr(colMajorDims);

    % coordinates (stored as netCDF type "variable")
    var_time = netcdf.defVar(ncid, 'time', 'float', dim_time);
    var_length = netcdf.defVar(ncid, 'length', 'float', dim_length);
    var_width = netcdf.defVar(ncid, 'width', 'float', dim_width);

    % store variables at the top level (add as needed)
    var_eta = netcdf.defVar(ncid, 'eta', 'float', rowMajorDims);
    var_discharge = netcdf.defVar(ncid, 'discharge', 'float', rowMajorDims);
    var_B = netcdf.defVar(ncid, 'channelwidth', 'float', rowMajorDims);
    var_H = netcdf.defVar(ncid, 'channeldepth', 'float', rowMajorDims);

    % store metadata variables as meta (add as needed)
    metaid = netcdf.defGrp(ncid, 'meta');
    var_H_SL = netcdf.defVar(metaid, 'H_SL', 'float', dim_time);
    var_dx = netcdf.defVar(metaid, 'dx', 'float', []);
    var_Qw = netcdf.defVar(metaid, 'Qw', 'float', []);
    var_Qs = netcdf.defVar(metaid, 'Qs', 'float', []);

    % now done creating things in the file
    netcdf.endDef(ncid);
    
    %% place variables in the fields created

    % place the constant variables into the file
    netcdf.putVar(ncid, var_length, fliplr(frst.grid.y(:,1)'));
    netcdf.putVar(ncid, var_width, frst.grid.x(1,:));
    netcdf.putVar(ncid, var_time, timesToSave);

    % loop through each timeslice, load and save variables as needed to a
    % temp arrays.
    temp_etas = zeros(length(outputFileList), ySize, xSize);
    temp_discharges = zeros(length(outputFileList), ySize, xSize);
    temp_Bs = zeros(length(outputFileList), ySize, xSize);
    temp_Hs = zeros(length(outputFileList), ySize, xSize);
    temp_H_SLs = zeros(1, length(outputFileList));
    for i=1:length(outputFileList)
        % load
        ith = load(fullfile(outputFileList(i).folder, outputFileList(i).name));

        % get variables
        i_eta = ith.grid.z;
        i_discharge = ith.grid.Qw;
        i_B = ith.grid.B;
        i_H = ith.grid.H;

        i_H_SL = ith.grid.oceanLevel;

        % do any variables need special treatment?
        %   hint: cannot store Inf in the netcdf file
        %   hint: best to replace NaN with 0, if possible
        i_B(~isfinite(i_B)) = 0;
        i_H(~isfinite(i_H)) = 0;

        % store each variable.
        %   NOTE: must permute each spatial variable,
        %   so that time is first dimension!
        temp_etas(i, :, :) = reshape(i_eta, 1, ySize, xSize);
        temp_discharges(i, :, :) = reshape(i_discharge, 1, ySize, xSize);
        temp_Bs(i, :, :) = reshape(i_B, 1, ySize, xSize);
        temp_Hs(i, :, :) = reshape(i_H, 1, ySize, xSize);

        temp_H_SLs(i) = i_H_SL;

    end

    % permute variables to expected ordering
    temp_etas  = permute(temp_etas, [3 2 1]);
    temp_discharges  = permute(temp_discharges, [3 2 1]);
    temp_Bs  = permute(temp_Bs, [3 2 1]);
    temp_Hs = permute(temp_Hs, [3 2 1]);

    % write out the eta field as time
    netcdf.putVar(ncid, var_time,timesToSave);
    netcdf.putVar(ncid, var_eta, temp_etas);
    netcdf.putVar(ncid, var_discharge, temp_discharges);
    netcdf.putVar(ncid, var_B, temp_Bs);
    netcdf.putVar(ncid, var_H, temp_Hs);

    netcdf.putVar(metaid, var_H_SL, temp_H_SLs);
    netcdf.putVar(metaid, var_dx, ith.grid.dx);
    netcdf.putVar(metaid, var_Qw, ith.grid.Qw(ith.grid.inletCell));
    netcdf.putVar(metaid, var_Qs, ith.grid.Qs_in(ith.grid.inletCell));

    % close the file
    netcdf.close(ncid)
 
end