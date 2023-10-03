function data = loadGamryDTA(filename, varargin)
%LOADGAMRYDTA Extract tables and metadata from a Gamry .DTA file.
%
% data = LOADGRAMRYDTA(filename) reads the contents of the Gamry DTA file
%   specified by FILENAME, extracts the metadata (tag, title, date, time),
%   and converts tabular data to MATLAB tables. DATA is a structure with
%   fields 'meta' and 'tables', structures containing the metadata and 
%   the tables (keyed by name) respectively.
%
% data = LOADGRAMRYDTA(filename,'NotesMode','KeyValue') also parses the
%   the notes section of the DTA file, which should take the following
%   key-value format:
%      KEY1: VALUE1
%      KEY2: VALUE2
%      ...
%   The DATA structure will contain an additional 'notes' field, a 
%   structure mapping each key to the associated value given in
%   the notes section of the DTA file. Numeric values are automatically
%   detected and converted to type double.
%
% -- Changelog --
% 2023.09.30 | Fix issue accessing table columns by name | Wes H.
% 2023.01.14 | Created | Wesley Hileman <whileman@uccs.edu>

p = inputParser;
p.addParameter('NotesMode','off', ...
    @(x)strcmpi(x,'off')||strcmpi(x,'KeyValue'));
p.parse(varargin{:});
notesMode = p.Results.NotesMode;

TAB = char(9);  % Horizontal tab character (used as delimiter in DTA files)
S_DEFAULT = 0;
S_TABLE = 1;
S_NOTES = 2;

S = readlines(filename);
data = struct;
state = S_DEFAULT;
nline = 1;
tabname = "";
tablen = 0;
noteslen = 0;

while nline <= length(S)
    if state == S_DEFAULT
        line = S(nline);
        parts = split(line,TAB);
        if strcmpi(parts(1),"TAG")
            data.meta.tag = parts(2);
        elseif strcmpi(parts(1),"TITLE")
            data.meta.title = parts(3);
        elseif strcmpi(parts(1),"DATE")
            data.meta.date = parts(3);
        elseif strcmpi(parts(1),"TIME")
            data.meta.time = parts(3);
        elseif strcmpi(parts(1),"PSTAT")
            data.meta.pstat = parts(3);
        elseif strcmpi(parts(1),"NOTES") && ~strcmpi(notesMode,"off")
            state = S_NOTES;
            noteslen = str2double(parts(3));
        elseif length(parts) > 1 && strcmpi(parts(2),"TABLE")
            state = S_TABLE;
            tabname = parts(1);
            % Length property not reliable, get length manually instead
            % (table lines are indented with one tab).
            tablen = 0;
            while true
                l = char(S(nline+tablen+1));
                if isempty(l) || l(1) ~= TAB
                    break;
                end
                tablen = tablen + 1;
            end
            if tablen < 2
                % Need at least two lines for header and units.
                error(['Malformed table: ' tabname]);
            end
            tablen = tablen - 2;  % Header and units take first two lines.
        end
        nline = nline + 1;
    elseif state == S_TABLE
        labels = split(strtrim(S(nline)),TAB)';
        units = split(strtrim(S(nline+1)),TAB)';
        tabdata = arrayfun( ...
            @(x)split(strtrim(x),TAB)', ...
            S(nline+2:nline+2+tablen-1), ...
            'UniformOutput',false);
        tabdata = vertcat(tabdata{:});
        tabdata = num2cell(tabdata,1);
        for col = 1:length(tabdata)
            asdouble = double(tabdata{col});
            if ~any(isnan(asdouble))
                % Use type double for this column.
                tabdata{col} = asdouble;
            end
        end
        tab = table(tabdata{:},'VariableNames',labels);
        tab.Properties.VariableUnits = units;

        % 2023.09.30 Not sure why this makes a difference, but it prevents 
        % errors sometimes when accessing columns by name. Happens at least
        % for large tables produced by GITT. MATLAB version R2021b.
        tab = tab(:,:);

        data.tables.(tabname) = tab;
        state = S_DEFAULT;
        tabname = "";
        tablen = 0;
        nline = nline + tablen + 2;
    elseif state == S_NOTES
        notesdata = arrayfun( ...
            @(x)split(strtrim(x),': ')', ...
            S(nline:nline+noteslen-1), ...
            'UniformOutput',false);
        for l = 1:length(notesdata)
            line = notesdata{l};
            if length(line)<2
                continue;
            end
            key = line(1);
            value = line(2);
            valueAsDouble = str2double(value);
            if ~isnan(valueAsDouble)
                value = valueAsDouble;
            end
            data.notes.(key) = value;
        end
        state = S_DEFAULT;
        noteslen = 0;
    end
end

end