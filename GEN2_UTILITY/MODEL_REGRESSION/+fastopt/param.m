function p = param(varargin)
    %PARAM Specify a parameter for optimization.
    %
    % p = PARAM(...) creates a specification for a parameter with
    %   fields set by name-value pairs in (...), similar to constructing a
    %   structure using the MATLAB built-in `struct` command.
    %
    % The available options are:
    %
    % - fix=VALUE: fix the value of the parameter to VALUE 
    %     (do not optimize)
    % - fixmask=BITMASK: for vector parameters, enables fixing individual
    %     elements of the vector. BITMASK is a logical vector of the same
    %     size as the value supplied for the `fix` parameter. Elements with
    %     false value are fixed to the corresponding element in the
    %     fix VALUE vector.
    % - len=LENGTH: for vector parameters, sets length as LENGTH 
    %     (default to length 1 scalars)
    % - logscale={true|false}: Default false. If true, store the parameter
    %     internally on a log10 scale (may help optimizers converge).
    % - tempfcn={'fix'|'lut'|'Eact'}: Default 'fix'.
    %     When 'fix', the parameter is constant with temperature
    %     When 'lut', the parameter takes on a different value for each
    %       temperature
    %     When 'Eact', the parameter follows an Ahhrenius relationship with
    %       temperature governed by an additional activation energy (Eact)
    %       parameter
    % - tempcoeff={'+'|'-'}. Default '+'. 
    %     When '+', use standard Ahhrenius relationship to temperature-
    %       correct the parameter (positive temperature coefficient). 
    %     When '-', use inverted Ahhrenius relationship to temperature-
    %       correct the parameter (negative temperature coefficient).

    % Create structure from suppied parameters.
    p = struct(varargin{:});

    % Set defaults / enforce integrity.
    if isfield(p,'fix'),         p.len = length(p.fix);         end
    if ~isfield(p,'len'),        p.len = 1;                     end
    if ~isfield(p,'fixmask'),    p.fixmask = false(p.len,1);    end 
    if ~isfield(p,'logscale'),   p.logscale = false;            end
    if ~isfield(p,'tempfcn'),    p.tempfcn = 'fix';             end
    if ~isfield(p,'tempcoeff'),  p.tempcoeff = '+';             end

    % Indicates this structure should not be flattened by FLATTENSTRUCT()
    p.noflatten__ = true;
end

