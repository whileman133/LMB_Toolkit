function builtstudy = build(studyname, currents, t0, tf, fs, varargin)
%BUILD Pre-process raw pulse-resistance data collected by a Gamry 
% potentiostat into a format usable by other scripts.
%
% INPUT:
% studyname = Path to directory containing the data folders for a SINGLE study.
% currents = vector of ordered pulse currents [A].
% t0 = time at onset of pulse [s].
% tf = time at offset of pulse [s].
% fs = nominal sampling rate [Hz].
%
% OPTIONAL KEYWORD ARGUMENTS:
% rawdir = directory in which raw study data resides (default PRROOT/raw)
% builtdir = directory in which built study data resides (default PRROOT/built)
%   This is the folder in which processed data is saved.
% only = cell array of the test names to include in the study (default all)
% verbose = logical, whether or not to show status messages in command
%   window (default true)
% outputname = name to give to the output study (default `studyname`)
% merge = boolean, when true merge with existing study (default false)
%
% studydir Directory Structure
%    Study1                     (Study or Cell folder - can correspond to a single cell, or multiple cells)
%    +---- Cell123_P25          (Cell or Temperature folder, P=+, N=- [Celcius], test name before underscore optional)
%    +-------- SOC_95-N5-35     (SOC folder, initial-increment-final, for increment P=+ and N=-)
%    +----------- ***.DTA       (DTA files output by the Gamry sequence)
%    +-------- SOC_30-N5-10
%    +----------- ***.DTA
%    +---- Cell456_N10
%    +-------- SOC_95-N5-35
%    +----------- ***.DTA  
%
% OUTPUT:
% study = `study` variable, which is a pr.BuiltStudy object with data fields. 
%   Access like a structure.

% Fetch root PR data directory.
PRROOT = TB.const.SDPROOT;

p = inputParser;
p.addOptional('only',{},@iscell);
p.addOptional('verbose',true,@islogical);
p.addOptional('rawdir',fullfile(PRROOT, 'labdata', 'raw'),@isfolder);
p.addOptional('builddir',fullfile(PRROOT, 'labdata', 'built'),@isfolder);
p.addOptional('outputname',studyname,@ischar);
p.addOptional('merge',false,@islogical);
p.parse(varargin{:});
whitelist = p.Results.only;
verbose = p.Results.verbose;
rawdir = p.Results.rawdir;
builddir = p.Results.builddir;
outputname = p.Results.outputname;
merge = p.Results.merge;

% Iterate cell/temperature folders.
tempfolders = getFolders(fullfile(rawdir,studyname));
tests = pr.Test.empty;
idxTest = 1;
for idxTemp = 1:length(tempfolders)
    tempfolder = tempfolders(idxTemp);
    temppath = fullfile(tempfolder.folder, tempfolder.name);

    % Determine temperature and test name from folder name.
    parts = split(tempfolder.name, '_');
    if length(parts) == 1
        temperature = getDouble(parts{1});
        testname = studyname;
    else
        testname = parts{1};
        temperature = getDouble(parts{2});
    end

    if verbose
        fprintf('Processing %s (T=%.0fC)... ', testname, temperature);
    end
    % Ensure this test is in the whitelist; skip if not.
    if ~isempty(whitelist) && ~any(strcmp(whitelist,testname))
        if verbose
            fprintf('SKIPPING, not in whitelist\n');
        end
        continue;
    end
    if verbose
        fprintf('\n');
    end

    % Iterate SOC folders - collect and sort pulse data grouped by SOC span.
    socfolders = getFolders(temppath);
    nSOC = length(socfolders);
    socGroups = struct( ...
        'sortindex',cell(1,nSOC),'soc',cell(1,nSOC),'path',cell(1,nSOC));
    for idxSOC = 1:nSOC
        socfolder = socfolders(idxSOC);
        socpath = fullfile(socfolder.folder, socfolder.name);

        % Determine SOC range from folder name.
        parts = split(socfolder.name, '_');
        socparts = split(parts{2}, '-');
        initial = str2double(socparts{1});
        increment = getDouble(socparts{2});
        final = str2double(socparts{3});

        socGroups(idxSOC).sortindex = initial;
        socGroups(idxSOC).soc = initial:increment:final;
        socGroups(idxSOC).path = socpath;
    end
    [~, idxSort] = sort([socGroups.sortindex],'descend');
    socGroups = socGroups(idxSort);

    % Create pulse objects.
    QtopAh = 0;
    pulses = [];
    for group = socGroups
        QdisAh = QtopAh + getDOD(group.path,group.soc);
        QtopAh = QdisAh(end);
        pulses = [pulses getPulses(group.path, group.soc, currents, QdisAh, t0, tf, fs, verbose)];
    end

    tests(idxTest) = pr.Test(testname, temperature, pulses);
    idxTest = idxTest + 1;
end % temperature

builtstudy = pr.BuiltStudy(outputname, tests);

% Determine save file name.
if ~isfolder(builddir), mkdir(builddir); end
studyFile = fullfile(builddir, sprintf('%s.mat', outputname));

% Merge studies if requested.
if (merge && isfile(studyFile))
    if verbose
        fprintf('\nLoading existing study %s for merge... ', studyFile);
    end
    studydata = load(studyFile,'builtstudy');
    oldstudy = studydata.builtstudy;
    if verbose
        fprintf('done!\n');
        fprintf('Merging studies... ');
    end
    builtstudy = pr.BuiltStudy.merge(outputname,builtstudy,oldstudy);
    if verbose
        fprintf('done!\n');
    end
end

% Save new study.
if verbose
    fprintf('\nSaving results to %s... ', studyFile);
end
save(studyFile, 'builtstudy', '-v7.3');
if verbose
    fprintf('done!\n');
    fprintf('\nTo fetch built study from disk, use:\n');
    fprintf('pr.load(''%s'',''built'')\n', ...
        builtstudy.name);
end

end

function pulses = getPulses(dtadir, socs, currents, QdisAh, t0, tf, fs, verbose)
    %GETPULSES Extract pulse data from DTA files.

    % Iterate DTA files containing pulse data.
    dtafiles = getFiles(dtadir);
    [~, idxSort] = sort({dtafiles.name});
    dtafiles = dtafiles(idxSort);
    pulses = pr.Pulse.empty; idxPulse = 1; 
    maxIdxIChg = 0; maxIdxIDis = 0;
    idxGroupIChg = 0; offsetIdxIChg = 0;
    idxGroupIDis = 0; offsetIdxIDis = 0;
    for idxFile = 1:length(dtafiles)
        dtafile = dtafiles(idxFile);
        dtapath = fullfile(dtafile.folder, dtafile.name);
        
        if dtafile.bytes == 0
            % Empty file, skip processing.
            continue;
        end

        match = regexp( ...
            dtafile.name, ...
            'CHRONOP(?<sign>DIS|CHG)(?<idxGroupI>[0-9]+)_#(?<idxSOC>[0-9]+)_#(?<idxI>[0-9]+)_#(?<idxRepeat>[0-9]+)', ...
            'names');
        if isempty(match)
            % Not a pulse file, skip processing.
            continue;
        end

        idxSOC = str2double(match.idxSOC);
        idxRepeat = str2double(match.idxRepeat);

        if idxSOC > length(socs)
            % Do not process - out of specified SOC range.
            continue;
        end

        % Resolve pulse-current index.
        newIdxGroupI = str2double(match.idxGroupI);
        if strcmp(match.sign, 'CHG')
            if newIdxGroupI > idxGroupIChg
                offsetIdxIChg = maxIdxIChg;
                idxGroupIChg = newIdxGroupI;
                maxIdxIChg = 0;
            elseif newIdxGroupI < idxGroupIChg
                error('Pulse-current group index should be ascending.');
            end
            idxI = offsetIdxIChg + str2double(match.idxI);
        else
            if newIdxGroupI > idxGroupIDis
                offsetIdxIDis = maxIdxIDis;
                idxGroupIDis = newIdxGroupI;
                maxIdxIDis = 0;
            elseif newIdxGroupI < idxGroupIDis
                error('Pulse-current group index should be ascending.');
            end
            idxI = offsetIdxIDis + str2double(match.idxI);
        end

        % Resolve SOC and current values.
        soc = socs(idxSOC);
        current = currents(idxI);
        if strcmp(match.sign, 'CHG')
            % Current sign is negative for charge using our convention.
            current = -current;
        end
        dod = QdisAh(idxSOC);

        % Extract pulse data.
        [t, i, v] = getTIVFromDTA(dtapath);
        i = -i;  % Gamry uses opposite sign convention for current!
        pulses(idxPulse) = ...
            pr.Pulse(t, v, i, soc, current, idxRepeat, dod, t0, tf, fs);

        % Status update.
        if verbose
            fprintf('%-30s%10d%10.3f%10d%10d of %10d\n', ...
                dtafile.name, soc, current, idxRepeat, idxFile, ...
                length(dtafiles));
        end

        % Store last pulse-current index for possible use next iteration.
        if strcmp(match.sign, 'CHG')
            maxIdxIChg = max(idxI, maxIdxIChg);
        else
            maxIdxIDis = max(idxI, maxIdxIDis);
        end
        
        idxPulse = idxPulse + 1;
    end % DTA files
end

function QdisAh = getDOD(dtadir, socs)
    %GETDOD Compute the total Ah discharged from the cell at each SOC
    %  setpoint from discharge DTA files.

    % Fetch discharge files.
    dtafiles = getFiles(dtadir);
    dtafiles = dtafiles(startsWith({dtafiles(:).name},'PWRDISCHARGE'));

    % Initial discharge before main test sequence (used to offset SOC
    % before test starts). Optional - may not exist.
    initial = strcmp({dtafiles(:).name},'PWRDISCHARGE_INITIAL.DTA');
    initialDis = dtafiles(initial);
    hasInitial = ~isempty(initialDis);

    % Main discharge steps executed during the test sequence.
    mainDis = dtafiles(~initial);

    QdisAh = zeros(length(socs),1);
    for idxFile = 1:length(mainDis)
        dtafile = mainDis(idxFile);
        dtapath = fullfile(dtafile.folder, dtafile.name);
        
        if dtafile.bytes == 0
            % Empty file, skip processing.
            continue;
        end

        % Determine SOC index from filename.
        parts = split(dtafile.name, '_'); 
        strIdx = parts{2};
        idxSOC = str2double(strIdx(2:end-4));

        if idxSOC > length(socs)
            % Do not process - out of specified SOC range.
            continue;
        end

        % Extract pulse data.
        [t, i, ~] = getTIVFromDTA(dtapath);
        i = -i;  % Gamry uses opposite sign convention for current!

        % Calculate total charge removed from the cell.
        QdisAh(idxSOC) = trapz(t,i)/3600;
    end

    % Make DOD cumulative.
    QdisAh = cumsum(QdisAh);

    % Apply initial offset if needed.
    if hasInitial
        [t, i, ~] = getTIVFromDTA(fullfile(initialDis.folder, initialDis.name));
        i = -i;  % Gamry uses opposite sign convention for current!
        QdisAh = QdisAh + trapz(t,i)/3600;
    end
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

function num = getDouble(string)
    %GETNUMBER Convert a string to number.
    if string(1)=='N'
        num=-str2double(string(2:end));
    else
        num=+str2double(string(2:end));
    end
end

function [time, current, voltage] = getTIVFromDTA(dtafilepath)
    %GETTIVFROMDTA Extract time, current, voltage data from a Gamry DTA file.

    % Readtable makes this relatively simple.
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