function builtstudy = build(studyname, Tref, Vmin, Vmax, varargin)
%BUILD Pre-process raw OCP data collected by a Gamry potentiostat 
% into a format usable by other scripts and save to disk.
%
% INPUT:
% studyname = Name of the study folder containing the raw data.
% Tref = Reference temperature for calibration scripts [C].
% Vmin, Vmax = Cell voltage limits [V].
%
% OPTIONAL KEYWORD ARGUMENTS:
% rawdir = directory in which raw study data resides (default OCPROOT/raw)
% builtdir = directory in which built study data resides (default OCPROOT/built)
%   This is the folder in which processed data is saved.
% only = cell array of the test names to include in the study (default all)
% verbose = logical, whether or not to show status messages in command
%   window (default true)
%
% celldir Directory Structure
%    Study1_0C01                     (Study or Cell folder - can correspond to a single cell, or multiple cells - C-rate after underscore - "C" occupies position of decimal point)
%    +---- Cell123_P25                (Cell or Temperature folder, P=+, N=- [Celcius], test name before underscore optional)
%    +-------- Script1                (OCP Script folder)
%    +------------Step 1.DAT          (Gamry data files to be concatenated,
%    +------------Step 2.DAT           sorted by name ascending)
%    +------------Step 3.DAT
%    +-------- Script2
%    +-------- Script3
%    +-------- Script4
%    +---- Cell456_N10
%    +-------- Script1
%    +-------- Script2
%    +-------- Script3
%    +-------- Script4
%    .
%    .
%    .
%
% OUTPUT:
% study = an ocp.Study object with data fields. Access like a structure. 
%   See ocp.Study.m, ocp.Test.m, and ocp.Script.m
%
% Also saves the study to disk.

% Fetch root OCP data directory.
OCPROOT = TB.const.OCPROOT;

p = inputParser;
p.addOptional('only',{},@iscell);
p.addOptional('verbose',true,@islogical);
p.addOptional('rawdir',fullfile(OCPROOT, 'labdata', 'raw'),@isfolder);
p.addOptional('builddir',fullfile(OCPROOT, 'labdata', 'built'),@isfolder);
p.parse(varargin{:});
whitelist = p.Results.only;
verbose = p.Results.verbose;
rawdir = p.Results.rawdir;
builddir = p.Results.builddir;

% Determine C-rate and study name from study folder name.
parts = split(studyname, '_');
if length(parts) == 1
    Crate = nan;
else
    Crate = getCrate(parts{2});
end

% Iterate temperature folders.
tempfolders = getFolders(fullfile(rawdir, studyname));
tests = ocp.Test.empty(length(tempfolders),0);
idxTest = 1;
for idxTemp = 1:length(tempfolders)
    tempfolder = tempfolders(idxTemp);
    temppath = fullfile(tempfolder.folder, tempfolder.name);

    % Determine temperature and test name from folder name.
    parts = split(tempfolder.name, '_');
    if length(parts) == 1
        temperature = getTemperature(parts{1});
        testname = studyname;
    else
        testname = parts{1};
        temperature = getTemperature(parts{2});
    end

    % Ensure this test is in the whitelist; skip if not.
    if verbose, fprintf('Processing %s (T=%.0fC)\n', testname, temperature); end
    if ~isempty(whitelist) && ~any(strcmp(whitelist,testname))
        if verbose, fprintf('\tSKIPPING, not in whitelist\n'); end
        continue;
    end

    % Iterate script folders.
    scriptfolders = getFolders(temppath);
    scripts = ocp.Script.empty(length(scriptfolders),0);
    for idxScript = 1:length(scriptfolders)
        scriptfolder = scriptfolders(idxScript);
        scriptpath = fullfile(scriptfolder.folder, scriptfolder.name);
        scriptname = lower(scriptfolder.name);
        idxscriptnbr = regexp(scriptname, '[0-9]+$');
        scriptnbr = str2double(scriptname(idxscriptnbr:end));

        if verbose, fprintf('\t%s', scriptname); end

        % ! Important. Gamry does not include the initial Time=0 sample.
        % (See also back interpolation below.)
        time = 0; 
        current = []; voltage = []; step = [];

        % Iterate script DTA files, concatenate data into contiguous
        % vectors.
        dtafiles = getFiles(scriptpath);
        dtafiles = sort({dtafiles.name});
        for idxFile = 1:length(dtafiles)
            dtafile = dtafiles{idxFile};
            dtapath = fullfile(scriptpath, dtafile);
            
            [t, i, v] = getTIVFromDTA(dtapath);

            % Sometimes Gamry provides one or more sample points with
            % the same time index (usu. at the end of the script). We
            % cannot allow duplicate data at the same time instant; 
            % remove them.
            [t, idxUniq] = unique(t,'first');
            i = i(idxUniq);
            v = v(idxUniq);

            time = [time; time(end)+t];
            current = [current; i];
            voltage = [voltage; v];
            step = [step; idxFile*ones(size(t))];

            if verbose, fprintf('.'); end
        end % DTA file

        % ! Important. Gamry does not include the initial Time=0 sample.
        % We back-interpolate the first recorded samples.
        current = [current(1); current];
        voltage = [voltage(1); voltage];
        step = [step(1); step];

        % Gamry uses opposite sign convention for current than ours!
        % Positive for anodic/oxidation current at the working
        % electrode, which corresponds to charging.
        current = -current;

        % Determine temperature at which script is performed.
        if any(strcmp(scriptname, {'script2', 'script4'}))
            % Calibration scripts performed at reference temperature.
            scriptT = Tref;
        else
            scriptT = temperature;
        end

        scripts(idxScript) = ocp.Script( ...
            scriptname, scriptnbr, scriptT, time, current, voltage, step);

        if verbose, fprintf(' OK\n'); end
    end % script

    tests(idxTest) = ocp.Test(testname, temperature, Vmin, Vmax, scripts);
    idxTest = idxTest + 1;
end % temperature

builtstudy = ocp.BuiltStudy(studyname, Tref, Crate, tests);

% Save processed data to disk.
if ~isfolder(builddir), mkdir(builddir); end
studyFile = fullfile(builddir, sprintf('%s.mat', builtstudy.name));
if verbose
    fprintf('\nSaving results to %s... ', studyFile);
end
save(studyFile, 'builtstudy');
if verbose
    fprintf('done!\n');
    fprintf('\nTo fetch built study from disk, use:\n');
    fprintf('ocp.load(''%s'',''built'')\n', ...
        builtstudy.name);
end

% Output results.
if verbose
    fprintf('\nSTUDY: %s\n', builtstudy.name);
    fprintf('%-20s%-10s%-10s%-10s\n', 'Test', 'T [degC]', 'QAh', 'eta');
    for test = builtstudy.tests
        fprintf('%-20s%-10.2f%-10.4f%-10.4f\n', ...
            test.name, test.temp, test.QAh, test.eta);
    end % tests
end

function folders = getFolders(path)
    %GETFOLDERS Get folders within a directory.
    nodes = dir(path);
    nodes = nodes([nodes.isdir]);
    nodes = nodes(~ismember({nodes(:).name},{'.','..'}));
    nodes = nodes(~startsWith({nodes(:).name},{'.'}));
    folders = nodes;
end

function files = getFiles(path)
    %GETFOLDERS Get files within a directory.
    nodes = dir(path);
    nodes = nodes(~[nodes.isdir]);
    nodes = nodes(~startsWith({nodes(:).name},{'.'}));
    files = nodes;
end

function Tnum = getTemperature(Tstr)
    %GETTEMPERATURE Convert a temperature string to number.
    Tnum = nan;
    if Tstr(1)=='N',Tnum=-str2double(Tstr(2:end));end
    if Tstr(1)=='P',Tnum=+str2double(Tstr(2:end));end
end

function Cnum = getCrate(Cstr)
    %GETCRATE Convert a C-rate string to number.
    cleaned = strrep(Cstr,'C','.');  % Re-insert decimal place.
    Cnum=str2double(cleaned);
end

function [time, current, voltage] = getTIVFromDTA(dtafilepath)
    %GETTIV Extract time, current, voltage data from a Gamry DTA file.

    % Readtable makes this easy!
    tab = readtable(dtafilepath, 'FileType', 'text');

    try
        time = tab.T;
    catch
        disp(tab)
        error('Problem reading data from %s', dtafilepath);
    end

    if any(strcmp('Im', tab.Properties.VariableNames))
        current = tab.Im;
    else
        % Assume this is an OCP measurement step; zero current flow.
        current = zeros(size(time));
    end
    voltage = tab.Vf;
end

end