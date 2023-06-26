function abspath = getAbsPath(relpath)
    % Locate this function.
    [fcndir,~,~] = fileparts(mfilename('fullpath'));
    % Locate the data directory relative to this function.
    abspath = fullfile(fcndir,'..',relpath);
end