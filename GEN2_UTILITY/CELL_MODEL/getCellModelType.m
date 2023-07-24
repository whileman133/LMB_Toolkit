function inModelType = getCellModelType(inModel,throwError)
%GETCELLMODELTYPE Determine the type of the input cell model.

if ~exist('throwError','var')
    throwError = true;
end

inModelType = [];
if isfield(inModel,'function')
    % Gen1 model: Legacy lumped-parameter model (LLPM).
    inModelType = 'LLPM';
elseif isfield(inModel,'type__')
    if strcmp(inModel.type__,'cellModel')
        % Gen2 model: P2DM, WORM, or RLWORM.
        inModelType = inModel.metadata.cell.type;
    end
end

if isempty(inModelType) && throwError
    error(['Unable to determine input model type.\n' ...
        'Should be Gen2 P2DM, WORM, or RLWORM; or Gen1 LLPM.']);
end

end