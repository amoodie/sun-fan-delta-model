function newAvulsions = avulsionCheck(grid,beta)
% avulsionCheck.m: Identifies new avulsion sites. For each channel cell, look up the local bed 
% elevation, channel depth, and downstream slope. Compare to elevation for
% neighboring cells and distance to those cells. Flag for avulsion
% if criterion is met. I imposed an additional constraint that a
% cell can only flow into 2 additional cells, so if the number of cells indicated grid.flowsTo is
% already 2 or more, no new avulsions can be made.

% initialize structure array to store row and column coordinates for new avulsion cells
newAvulsions.rSource = []; 
newAvulsions.cSource = [];
newAvulsions.indSource = [];
newAvulsions.rNew = [];
newAvulsions.cNew = [];
newAvulsions.indNew = [];
        
    % Check for new avulsions at each channel cell
    for k=1:numel(grid.channelFlag)
        if grid.channelFlag(k) && numel(grid.flowsTo{k})<2 % i.e., if it's a channel cell and flows to no more than 1 cell, then eligible for a new avulsion
            z_i = grid.z(k);
            H_ij = grid.H(k); % listed as "H_ij" in the paper, but I don't understand that notation as it seems to imply depth from i to j; easier to conceive as depth at i

            if isinf(H_ij)|| isnan(H_ij) || H_ij <= 0 
                % because depth is calculated using slope, and slope can be negative, the value of 
                % depth can be non-physical (Inf / NaN / negative). In that
                % case, it cannot be used in the avulsion criterion in
                % equation 13; I think a reasonable workaround is to
                % set H_ij to zero for this case so that flow depth
                % does not impact avulsion susceptibilty for this case.
                H_ij = 0;
            end

            S_ij = grid.S.alongFlow(k);
            % search nearest neighbor cells for potential avulsion
            % destinations. Skip any current channel cells or ocean
            % cells.
            [currentRow,currentCol] = ind2sub(grid.size,k);

            % search the 8 nearest neighbor cells
            rowSearch = currentRow + [-1 -1 -1 0 1 1 1 0]; % clockwise from NW
            colSearch = currentCol + [-1 0 1 1 1 0 -1 -1];
            L_ik = grid.dx*[sqrt(2) 1 sqrt(2) 1 sqrt(2) 1 sqrt(2) 1]; % distance from starting cell to search cell
            avulsionSusceptibilityIndex = nan(1,8); % NaN rather than 0 so that "0" doesn't accidentally become a maximum
            for l=1:numel(rowSearch)
                if rowSearch(l) < 1 || rowSearch(l) > grid.size(1) || colSearch(l)<1 || colSearch(l) > grid.size(2) || grid.channelFlag(rowSearch(l),colSearch(l))
                    continue % if search cell is already a channel cell or is off the grid then cannot avulse to there. Therefore, leave avulsion susceptibility index as NaN for this cell and continue to next loop iteration.
                else
                    z_k = grid.z(rowSearch(l),colSearch(l));
                    avulsionSusceptibilityIndex(l) = ((z_i-beta*H_ij) - z_k)/L_ik(l); % equation 13 LHS
                end
            end

            if any(avulsionSusceptibilityIndex > S_ij) % equation 13
                % select the neighboring cell with the greatest avulsion
                % susceptibility 
                [~,neighborAvulsionSelect] = max(avulsionSusceptibilityIndex);
                newAvulsions.rNew = [newAvulsions.rNew; rowSearch(neighborAvulsionSelect)];
                newAvulsions.cNew = [newAvulsions.cNew; colSearch(neighborAvulsionSelect)];
                indNew = sub2ind(grid.size,rowSearch(neighborAvulsionSelect),colSearch(neighborAvulsionSelect));
                newAvulsions.indNew = [newAvulsions.indNew; indNew];
                newAvulsions.rSource = [newAvulsions.rSource; currentRow];
                newAvulsions.cSource = [newAvulsions.cSource; currentCol];
                newAvulsions.indSource = [newAvulsions.indSource; k];
            end
        end
    end
end