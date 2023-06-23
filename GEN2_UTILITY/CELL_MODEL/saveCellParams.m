function saveCellParams(model, filename, varargin)
    %SAVECELLPARAMS Write a cell model to an Excel workbook.
    %
    %  SAVECELLPARAMS(model,filename) saves the cell model specified by
    %    MODEL to the Excel workbook located in the filesystem at FILENAME.
    %    Re-load the model from the filesystem using LOADCELLPARAMS.
    %
    %  SAVECELLPARAMS(...,'MSMR',false) saves the cell model without the
    %    MSMR parameters (if the model is an MSMR model), writing OCP and OCV
    %    tables instead. The MSMR parameters are discarded. 
    %    For non-MSMR models, there is no change.
    %
    %  SAVECELLPARAMS(...,'verbose',true) prints warning messages to the
    %    console indicating missing or extra fields in the template file.
    %    By default, these messages are suppressed.
    %
    % -- Limitations --
    % Make sure to conform to the following if you edit the cell model
    % (the inline function parameters) manually.
    %
    %   * Lookup Tables (LUTs): Use INTERP1 to perform the interpolation.
    %
    %   * Custom functions: Reference a function parameter in the function 
    %     body at least once. For example, the following will be recognized 
    %     as a custom function:
    %
    %        @(x,T)(rand()+0*x)
    %
    %     while this will not:
    %
    %        @(x,T)(rand())
    %
    %   * Constant Parameters: All parameters not recognized as either LUTs
    %     or custom functions are interpreted as constants. Given an inline
    %     function denoted fcn(x,T), the value of the constant is taken as
    %     fcn(0,0) (that is, the function is evalulated for x=0 and T=0 to
    %     determine the constant). The constant may be a vector or a
    %     scalar.
    %
    %   * Activation Energy: If a parameter has activation energy Ea,
    %     specify the temperature depedence *exactly* as follows to ensure
    %     Ea is reversed correctly:
    %
    %        <function>.*exp(Ea*(1/298.15-1/T)/8.31446)
    %
    %     where <function> is the inline function without the temperature
    %     dependence. Ea may be a scalar or a vector. LUTs, custom
    %     functions, and constants may all be associated with activation
    %     energy.
    %
    % -- Changelog --
    % 05.02.2022 | Wesley Hileman | Created

    SHEET = 'Parameters';  % Template sheet in which to store parameters.
    COLVALUE = 'C';        % Column in which to place parameter values.
    COLEACT = 'D';         % Column in which to place parameter activation energy (should be adjacent to COLVALUE).

    % Parse optional arguments.
    parser = inputParser;
    addParameter(parser,'MSMR',true,@islogical);
    addParameter(parser,'verbose',false,@islogical);
    parse(parser,varargin{:});

    % Decide which cell template to use.
    useMSMR = parser.Results.MSMR && model.MSMR;
    [fcndir,~,~] = fileparts(mfilename('fullpath'));
    if useMSMR
        template = fullfile(fcndir,'..','XLSX_CELLDEFS','cellTEMPLATE-MSMR.xlsx');
    else
        template = fullfile(fcndir,'..','XLSX_CELLDEFS','cellTEMPLATE.xlsx');
    end

    % Console output selector.
    verbose = parser.Results.verbose;

    % Prepare new workbook, determine parameter row positions in workbook.
    [succ,errmsg,errid] = copyfile(template,filename);
    if ~succ
        error('Error creating "%s": %s - %s',filename,errid,errmsg);
    end
    dict = buildParamDict(filename,SHEET);

    % Write general section.
    writeCells(model.name,filename,SHEET,dict.gen.name,COLVALUE);
    writeCells(1,filename,SHEET,dict.gen.lumped,COLVALUE);  % ! We always write the lumped parameters out.
    writeCells(double(useMSMR),filename,SHEET,dict.gen.MSMR,COLVALUE);

    % Write remaining sections.
    luts = struct('sheetname',{},'value1',{},'value2',{});  % Storage for LUTs to be written later.
    fcns = model.function;
    secs = fieldnames(fcns);
    for idxSec = 1:length(secs)
        sec = secs{idxSec};
        if ~isfield(dict,sec)
            if verbose
                warning(['Template is missing section "%s". ' ...
                         'Omitting from export.'],sec);
            end
            continue;
        end

        % Ensure template contains all parameters in model, generate
        % warnings if not.
        paramsmodel = fieldnames(fcns.(sec));
        paramstempl = fieldnames(dict.(sec));
        paramsmissing = setdiff(paramsmodel,paramstempl);
        paramsextra = setdiff(paramstempl,paramsmodel);
        if verbose
            if ~isempty(paramsmissing)
                warning(['Template is missing parameters "%s" in section ' ...
                         '"%s". Omitting from export.'],strjoin(paramsmissing,', '),sec);
            end
            if ~isempty(paramsextra)
                warning(['Template has extra parameters "%s" in section ' ...
                         '"%s".'],strjoin(paramsextra,', '),sec);
            end
        end

        % Build cell array to insert into section.
        % First column is parameter value, second column is Ea.
        % Note that the fields in MATLAB structures are ordered, so those
        % in the `dict` variable appear in the same order as they do in the
        % spreadsheet.
        firstrow = dict.(sec).(paramstempl{1}); 
        lastrow = dict.(sec).(paramstempl{end});
        seclength = lastrow - firstrow + 1;
        secdata = cell(seclength,2);
        for idxParam = 1:length(paramstempl)
            param = paramstempl{idxParam};
            if ~isfield(fcns.(sec),param)
                % Skip extra template parameters.
                continue;
            end
            idxRow = dict.(sec).(param) - firstrow + 1;

            % Convert function back to object.
            obj = reverseFcn(fcns.(sec).(param));
            obj.Ea = obj.Ea/1000; % Convert from J/mol to kJ/mol.

            switch(obj.type)
                case 'LUT'
                    % Save LUT for later write.
                    secname = [upper(sec(1)) sec(2:end)];
                    sheetname = [param secname];
                    luts(end+1) = ... 
                        struct('sheetname',sheetname,'value1',obj.value{1},'value2',obj.value{2});
                    % Create link to LUT.
                    secdata{idxRow,1} = ['#' sheetname];
                    secdata{idxRow,2} = obj.Ea;
                case 'function'
                    % Convert function definition to string
                    fcnStr = func2str(obj.value);
                    secdata{idxRow,1} = fcnStr;
                    secdata{idxRow,2} = obj.Ea;
                case 'vector'
                    if isscalar(obj.Ea)
                        Ea = obj.Ea;
                    else
                        Ea = mat2str(obj.Ea(:)');
                    end
                    secdata{idxRow,1} = mat2str(obj.value(:)');
                    secdata{idxRow,2} = Ea;
                case 'scalar'
                    secdata{idxRow,1} = obj.value;
                    secdata{idxRow,2} = obj.Ea;
            end % switch
        end % parameter loop

        % Write section to workbook.
        tab = cell2table(secdata);
        writeCells(tab,filename,SHEET,[firstrow lastrow],[COLVALUE COLEACT]);
    end % section loop

    % Write each LUT to workbook.
    if ~isempty(luts)
        for lut = luts
            tab = array2table([lut.value1 lut.value2]);
            writeCells(tab,filename,lut.sheetname,[1 length(lut.value1)],'AB');
        end
    end
end

function writeCells(tab, filename, sheet, row, col)
    %WRITECELLS Write data to a block of cells in an Excel workbook.

    if ~istable(tab)
        if ischar(tab); tab = {tab}; end
        tab = table(tab);
    end
    if isscalar(row) && isscalar(col)
        cells = sprintf('%s%d',col,row);
    elseif isscalar(row)
        cells = sprintf('%s%d:%s%d',col(1),row,col(2),row);
    elseif isscalar(col)
        cells = sprintf('%s%d:%s%d',col,row(1),col,row(2));
    else
        cells = sprintf('%s%d:%s%d',col(1),row(1),col(2),row(2));
    end

    try
      writetable(tab,filename, ...
          'WriteVariableNames',0, ...
          'Sheet',sheet, ...
          'Range',cells, ...
          'AutoFitWidth',false, ...
          'PreserveFormat',true);
    catch ME
      if ~strcmp(ME.identifier,'MATLAB:table:parseArgs:BadParamName')
        rethrow(ME);
      end
      % Older versions of MATLAB do not support preserving spreadsheet
      % formatting.
      warning(['Spreadsheet formatting will change. ' ...
               'Use MATLAB >= R2020a to preserve formatting.']);
      writetable(tab,filename, ...
        'WriteVariableNames',0, ...
        'Sheet',sheet, ...
        'Range',cells);
    end
end

function dict = buildParamDict(filename, sheet)
    %BUILDPARAMDICT Associate parameter names with row numbers in an Excel
    % sheet.
    %
    % dict = BUILDPARAMDICT(filename,sheet) associates the codenames in
    %  each section of the parameter table located in SHEET of the Excel 
    %  workbook at FILENAME with linenumbers. DICT is a structure
    %  containing sub-structures for each section. Each section structure
    %  maps parameter codenames to line numbers in the excel worksheet.

    % Number of lines to skip at start of each section 
    % (hidden line and table header).
    HEADSKIP = 3;

    [ind,data] = loadParamTable(filename,sheet);
    sections = fieldnames(ind);
    for idxSec = 1:length(sections)
        sec = sections{idxSec}; indsec = ind.(sec);
        secdata = data(indsec,:); 
        secdata = secdata(HEADSKIP+1:end,:);
        codenames = secdata(:,2);
        for k = 1:length(codenames)
            name = codenames{k};
            if ~ischar(name); continue; end  % Skip empty lines.
            dict.(sec).(name) = indsec(1) + HEADSKIP + k - 1;
        end
    end
end
