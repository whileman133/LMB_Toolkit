function addPaths()
    %ADDPATHS Add all LMB toolbox directories to the MATLAB path.
    addpath(TB.getAbsPath('UTILITY'));
    addpath(TB.getAbsPath('TFS'));
    addpath(TB.getAbsPath('XLSX_CELLDEFS'));
    addpath(TB.getAbsPath('MAT_CELLDEFS'));
end