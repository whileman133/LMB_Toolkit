function cellData = tfCommon(s,cellData)
% TFCOMMON
%
% -- Changelog --
% 2023.10.31 | Completely refactor for gen2 toolkit | Wes H.
% 2023.09.17 | Use tfLMB() to evalulate parameter values | Wes H.
% 2023.09.14 | Use Warburg parameters | Wesley Hileman <whileman@uccs.edu>

% first, check to see if we have already computed the relevant data...
% don't recompute unless necessary!
if ( ...
    isfield(cellData,'common') && ...
    isfield(cellData.common,'s') && ...
    isequal(cellData.common.s(:),s(:)) ...
)
    return;
end

% Fetch cell parameters from tfLMB function. This keeps all of the
% parameter evalulation in a single place. Modify tfLMB to change how
% the parameter values are evalulated.
tfData = tfLMB(s,cellData, ...
    'Calc11',true,'Calc22',false, ...
    'TdegC',cellData.const.T-273.15,'socPct',cellData.const.soc*100);
cellData.common.s = s;
cellData.common.tfData = tfData;

end % function