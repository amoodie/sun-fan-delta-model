function [avulsionCellInfo] = avulsionCheck(grid,beta,Qs_threshold,branchLimit)
% avulsionCheck.m: Identifies new avulsion sites. For each channel cell,
% look up the local bed elevation, channel depth, and downstream slope.
% Compare to elevation for neighboring cells and distance to those cells.
% Flag for avulsion if criterion is met. An additional constraint is
% that a cell can only flow into 2 additional cells, so if the number of
% cells indicated grid.flowsTo is already 2 or more, no new avulsions can
% be made. This is different from the Sun et al. (2002) model.

% initialize an array to store the cell indices where avulsion should happen
avulsionCellInfo = [];

    % determine all the channel locations to loop through
    channelInds = find(grid.channelFlag);

    L_ik = grid.dx*[sqrt(2); 1; sqrt(2); 1; sqrt(2); 1; sqrt(2); 1]; % distance from starting cell to search cell

    % Check for new avulsions at each channel cell
    for kk=1:numel(channelInds)

        k = channelInds(kk);

        if grid.flowsToCount(k) < branchLimit % i.e., if it flows to no more than branchLimit cells, then eligible for a new avulsion
            % get flow depth at current cell
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

            % get current down flow slope
            S_ij = grid.S.alongFlow(k);

            % determine avulsion susceptibility for ngbrs
            avulsionSusceptibilityIndex = NaN(8,1);
            %   only allowed to avulse from cells with sediment flux
            if grid.Qs_in(k) >= Qs_threshold
                % get idxs of neighbors
                nghbrs = grid.nghbrs(:, k);
                nghbrsValid = (nghbrs ~= 0);

                % compute avulsion susceptibilty
                avulsionSusceptibilityIndex(nghbrsValid) = ((grid.z(k) - (beta.*H_ij)) - grid.z(nghbrs(nghbrsValid))) ./ L_ik(nghbrsValid);
            end

            % check if cell is susceptible to avulsion!
            if any(avulsionSusceptibilityIndex > S_ij) % equation 13
                % which cell should the avulsion go into (most likely)
                %   we pass the likely info to avulsion router
                [~, nghbrToAvulseTo] = max(avulsionSusceptibilityIndex);
                avulsionInfo = [k, nghbrToAvulseTo];

                % flag this cell for an avulsion
                avulsionCellInfo = [avulsionCellInfo; avulsionInfo];
            end
        end
    end
end
