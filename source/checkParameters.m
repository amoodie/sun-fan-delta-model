function checkParameters(parameters)
% checkParameters.m: Checks that each parameter has the expected data type
% and if not, throws an error.

msg = 'Incorrect data type for input parameter';

% Check data type for each parameter
if ~ischar(parameters.runName), error(msg); end
if ~islogical(parameters.clobber), error(msg); end
if ~isnumeric(parameters.alpha_so),error(msg); end
if ~isnumeric(parameters.alpha_sa),error(msg); end
if ~isnumeric(parameters.alpha_r),error(msg); end
if ~isnumeric(parameters.alpha_b),error(msg); end
if ~isnumeric(parameters.tauStar_c),error(msg); end
if ~isnumeric(parameters.n),error(msg); end
if ~isnumeric(parameters.p),error(msg); end
if ~isnumeric(parameters.R),error(msg); end
if ~isnumeric(parameters.g),error(msg); end
if ~isnumeric(parameters.gamma),error(msg); end 
if ~isnumeric(parameters.lambda),error(msg); end
if ~isnumeric(parameters.beta),error(msg); end 
if ~isnumeric(parameters.Qw_inlet),error(msg); end
if ~isnumeric(parameters.Qs_inlet),error(msg); end
if ~isnumeric(parameters.Qw_threshold),error(msg); end
if ~isnumeric(parameters.Qw_mismatch_tolerance),error(msg); end
if ~isnumeric(parameters.D),error(msg); end
if ~isstruct(parameters.oceanLevel),error(msg); end
if ne(numel(parameters.oceanLevel.timeStart_yr),numel(parameters.oceanLevel.z))
    error('In oceanLevel, fields timeStart_yr and z must have the same number of elements');
end
if ~isstruct(parameters.grid), error(msg); end
if ~isnumeric(parameters.t), error(msg); end
if ~isnumeric(parameters.tStep_sec), error(msg); end
if ~isnumeric(parameters.tMax_yr), error(msg); end
if ~isnumeric(parameters.tMax_sec), error(msg); end
if ~isnumeric(parameters.startingTime), error(msg); end
if ~isnumeric(parameters.tSaveInterval_yr), error(msg); end
if ~isnumeric(parameters.secondsPerYear), error(msg); end
if ~isnumeric(parameters.tElapsedSinceSave_yr), error(msg); end
if ~isstruct(parameters.inlet), error(msg); end
if ~ischar(parameters.boundaryCondition), error(msg); end
if ~islogical(parameters.debugFigure), error(msg); end
if ~islogical(parameters.loadCheckpoint) && ~ischar(parameters.loadCheckpoint), error(msg); end
if ~ischar(parameters.outputDir), error(msg); end

end