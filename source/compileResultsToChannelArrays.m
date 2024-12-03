function compileResultsToChannelArrays(resultsFolder)

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

    % set up the output arrays
    flowsTo = zeros(ySize, xSize, 2);
    B = zeros(ySize, xSize);
    
    % which indexs do we want to plot?
    outputFileList = outputFileList(end);
    
    % loop through each timeslice, load and save to a temp array
    for i=1:length(outputFileList)
        % load
        ith = load(fullfile(outputFileList(i).folder, outputFileList(i).name));
       
        % grab the B array
        B(:) = ith.grid.B(:);
        
        % loop through every cell
        for j=1:xSize
            for k=1:ySize
                jkFlowsTo = ith.grid.flowsTo{j, k};
                
                if numel(jkFlowsTo) == 0
                    continue
                elseif numel(jkFlowsTo) == 1
                    flowsTo(j, k, 1) = jkFlowsTo;
                elseif numel(jkFlowsTo) == 2
                    flowsTo(j, k, :) = jkFlowsTo(:);
                else
                    error('bad formatting.')
                end
            end
        end
        
        % write out the mat file to load in python
        save(fullfile(outputFileList(i).folder, ['channelNetwork_','postdoc_sloped_nowater_short_time_50.00_yr.mat']), 'flowsTo', 'B')
        
    end

    
    
end