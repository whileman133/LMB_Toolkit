function inModelType = getCellModelType(inModel)
%GETCELLMODELTYPE Determine the type of the input cell model.

if isfield(inModel,'metadata') && ...
   isfield(inModel,'type__') && ... 
   strcmp(inModel.type__,'cellModel')
    % Gen2 model: P2DM, WORM, or RLWORM.
    inModelType = inModel.metadata.cell.type;
elseif isfield(inModel,'function')
    % Gen1 model: Legacy lumped-parameter model (LLPM).
    inModelType = 'LLPM';
else
    error(['Unable to determine input model type.\n' ...
        'Should be P2DM, WORM, RLWORM, or Gen1 LLPM.']);
end % if

end