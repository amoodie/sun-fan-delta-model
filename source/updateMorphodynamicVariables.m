function grid = updateMorphodynamicVariables(grid,alpha_b,alpha_r,alpha_sa,alpha_so,R,g,D,tauStar_c,n,p)
% updateMorphodyamicVariables.m: for each cell on the grid, updates channel
% width (grid.B), channel depth (grid.H), and sediment flux out of the cell
% (grid.Qs_out). 

    % calculate channel width (eqn. 9a), depth (eqn. 9b), sediment discharge
    % (eqn. 9c). Each of these is written in dimensionless form in the
    % corresponding equations, but most of their results plots (Fig. 2 to 7)
    % use dimensions. So here I calculate width, depth, and sediment discharge
    % in their dimensioned forms. Each quantity is calculated
    % simulataneously for all cells on the grid.

    % slope is defined as positive downhill.
    grid.B = (alpha_b^(-(3+2*p)/2)).*(alpha_r^-1).*(grid.S.alongFlow.^(1+p)).*((grid.Qw/(sqrt(g*D)*D^2))).*D; % width
    grid.H = (alpha_b./grid.S.alongFlow)*D; % depth
    % Note that given these dependencies on slope, both depth and width
    % can come out with non-physical values (i.e., negative or zero).
    % Those cases would also result in Qs_out = 0 based on the code
    % immediately below. The only other effect is that depth is used in
    % the avulsion check;  I think (but haven't double-check) cells with negative depth sd
    % cannot meet the avulsion criteria.

    grid.Qs_out = alpha_so*alpha_sa*(alpha_b^(-(3+2*p)/2))*(alpha_r^-1)*((alpha_b/R - tauStar_c)^n)*(grid.S.alongFlow.^(1+p)).*(grid.Qw/(sqrt(g*D)*D^2))*(sqrt(R*g*D)*D^2); % sediment discharge

    % set any points with Qs_out < 0 to Qs_out = 0 (can happen for
    % negative slope values (i.e., adverse slope)). In those cases,
    % no sediment should leave the cell.
    grid.Qs_out(grid.Qs_out<0) = 0;

    % Enforce no sediment flux out for any cells that do not flow to
    % other cells
    grid.Qs_out(cellfun(@isempty,grid.flowsTo)) = 0;

    % Enforce no sediment flux out the for any cells below sea
    % level. (This is likely redundant, as no ocean cells should flow
    % to any other cells). 
    grid.Qs_out(grid.oceanFlag) = 0;
        
    if any(isnan(grid.Qs_out(grid.channelFlag)))
        error('Unexpected NaN value in sediment flux calculation');
    end
end