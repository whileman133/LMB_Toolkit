function load(studyname, kind, varargin)
    %LOAD Load OCP data from disk.
    %
    % LOAD(studyname,'built') reads the built (i.e. preprocessed) 
    %   OCP data part of the study STUDYNAME from the default 
    %   directory and assigns it in the calling workspace. The data is not 
    %   re-loaded should it already exist in the calling workspace.
    %
    % LOAD(studyname,'fit') reads the fit (i.e. regressed) OCP
    %   models part of the study STUDYNAME from the default OCP "fit" data
    %   directory and assigns it in the calling workspace. 
    %   The data is not re-loaded should it already exist in the calling workspace.
    %
    % LOAD(...,'studydir',DIRECTORY) loads the built study data from the
    %   specified directory instead of the default location.
    %
    % LOAD(...,'verbose',false) suppresses all command-window output.
    %
    % LOAD(...,'varname',NAME) assigns the OCP data to the variable named
    %   NAME in the calling workspace instead of the default.

    parser = inputParser;
    addOptional(parser,'studydir','',@isfolder);
    addOptional(parser,'verbose',true,@islogical);
    addOptional(parser,'varname','',@ischar);
    parse(parser,varargin{:});

    if strcmp(kind, 'built')
        ocptype = 'ocp.BuiltStudy';
        defvarname = 'builtstudy';
    elseif strcmp(kind, 'fit')
        ocptype = 'ocp.FitStudy';
        defvarname = 'fitstudy';
    else
        error(['Unrecognized OCP data specifier. ' ...
            'Should be ''built'' or ''fit''.'])
    end
    
    % Determine directory in which processed study data resides.
    studydir = parser.Results.studydir;
    if isempty(studydir)
        studydir = fullfile(TB.const.OCPROOT,'labdata',kind);
    end

    % Set output level.
    verbose = parser.Results.verbose;

    % Determine variable name.
    varname = parser.Results.varname;
    if isempty(varname)
        varname = defvarname;
    end

    % Try to find existing data variable in the calling workspace.
    try
        data = evalin('caller',varname);
    catch
    end
    
    % Load study data if does not exist in calling workspace.
    if ~exist('data','var') || ~isa(data,ocptype) || ~strcmp(data.name,studyname)
        if verbose
            fprintf(['Loading %s OCP data for ''%s'' study. ' ...
                'This may take a few moments... '], kind, studyname);
        end
        studyFile = fullfile(studydir, sprintf('%s.%s', studyname, 'mat'));
        data = load(studyFile);
        fields = fieldnames(data);
        data = data.(fields{1});
        if verbose
            fprintf('done!\n');
        end
    elseif verbose
        fprintf( ...
            'Using previously-loaded %s data for ''%s'' study.\n', ...
            kind, studyname);
    end

    % Assign study in calling workspace.
    assignin('caller',varname,data);
end