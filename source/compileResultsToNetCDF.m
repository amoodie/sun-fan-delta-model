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

    dim_length = netcdf.defDim(ncid,'length',xSize);
    dim_width = netcdf.defDim(ncid,'width',ySize);
    dim_time = netcdf.defDim(ncid,'total_time',length(outputFileList));
         
    colMajorDims = [dim_time, dim_length, dim_width];
    rowMajorDims = fliplr(colMajorDims);
    
    var_x = netcdf.defVar(ncid,'x','float',[dim_length, dim_width]);
    var_y = netcdf.defVar(ncid,'y','float',[dim_length, dim_width]);
    var_time = netcdf.defVar(ncid,'time','float',dim_time);

    metaid = netcdf.defGrp(ncid,'meta');
    var_eta = netcdf.defVar(ncid,'eta','float',rowMajorDims);
    var_H_SL = netcdf.defVar(metaid,'H_SL','float',dim_time);
    
    netcdf.endDef(ncid);
    
    %% place variables in the fields created

    % place the constant variables into the file
    netcdf.putVar(ncid, var_x, permute(frst.grid.x/frst.grid.dx, [2, 1]));
    netcdf.putVar(ncid, var_y, permute(frst.grid.y/frst.grid.dx, [2, 1]));
    netcdf.putVar(ncid,var_time,timesToSave);
         
    % loop through each timeslice, load and save to a temp array (ts)
    ts = zeros(length(outputFileList), xSize, ySize);
    H_SL = zeros(1, length(outputFileList));
    for i=1:length(outputFileList)
        % load
        ith = load(fullfile(outputFileList(i).folder, outputFileList(i).name));
        
        zi = ith.grid.z;
        if i == 1
            H_SL(i) = 80;
        else
            H_SL(i) = ith.grid.oceanLevel;
        end
        
        % permute so time is first dimensions
        permuted = reshape(zi, 1, xSize, ySize);
        
        ts(i, :, :) = permuted;
    end
    
    % permute ts to rowMajorDims
    ts  = permute(ts, [3 2 1]);
    
    % write out the eta field as time
    netcdf.putVar(ncid, var_time,timesToSave);
    netcdf.putVar(ncid, var_eta, ts);
    netcdf.putVar(metaid, var_H_SL, H_SL);
    
    % close the file
    netcdf.close(ncid)
 
end