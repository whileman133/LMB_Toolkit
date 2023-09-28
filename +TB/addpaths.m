function addpaths(varargin)
    %ADDPATHS Add all LMB toolbox directories to the MATLAB path.

    parser = inputParser;
    parser.addOptional('prefix','gen2',@ischar);
    parser.parse(varargin{:});
    prefix = parser.Results.prefix;
    if ~isempty(prefix)
        prefix = [upper(prefix) '_'];
    end

    addpath(genpath(TB.getAbsPath([prefix 'UTILITY'])));
    addpath(genpath(TB.getAbsPath([prefix 'TFS'])));
    addpath(genpath(TB.getAbsPath([prefix 'XLSX_CELLDEFS'])));
    addpath(genpath(TB.getAbsPath([prefix 'MAT_CELLDEFS'])));
    addpath(genpath(TB.getAbsPath([prefix 'XLSX_INPUTS'])));
    addpath(genpath(TB.getAbsPath([prefix 'XLSX_XRACTRL'])));
end