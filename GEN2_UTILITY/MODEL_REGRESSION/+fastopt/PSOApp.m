classdef PSOApp < handle
%PSOAPP

properties(SetAccess=protected)
    arg           % structure of validated input arguments
    freeParam     % structure of parameters allowed to optimize over
    dof           % number of degrees of freedom in model
    mult          % multiplicity of the model (number of temperatures)
    TdegC         % temperature vector (degC)
    lb0           % lower bounds for reset
    ub0           % upper bounds for reset
    init0         % initial values for reset

    % UI components.
    fig           % figure window
    panel         % structure of panels of the figure window
                  % (.control, .parameter, .plot)
    btn           % structure of uibutton
    edt           % structure of uieditfield
    lab           % structure of uilabel
    sld           % structure of uislider
    status        % program status output (uilabel)

    % State variables.
    startFlag
    values
    lb
    ub
end

methods
    function obj = PSOApp(modelspec,init,lb,ub,varargin)
        parser = inputParser;
        parser.addRequired('modelspec',@isstruct);
        parser.addRequired('init',@isstruct);
        parser.addRequired('lb',@isstruct);
        parser.addRequired('ub',@isstruct);
        parser.addParameter('SwarmSize',500);
        parser.addParameter('SwarmMaximumIterations',100);
        parser.addParameter('FminMaximumIterations',1000);
        parser.addParameter('FminMaximumFunctionEvals',10000);
        parser.addParameter('PlotInitializeFcn', ...
            @(varargin)[],@(x)isa(x,'function_handle'));
        parser.addParameter('PlotUpdateFcn', ...
            @(varargin)[],@(x)isa(x,'function_handle'));
        parser.addParameter('WindowName','Particle Swarm Optimization');
        parser.parse(modelspec,init,lb,ub,varargin{:});
        arg = parser.Results;  % structure of validated arguments
        obj.arg = arg;
        obj.init0 = fastopt.flattenstruct(arg.init);
        obj.lb0 = fastopt.flattenstruct(arg.lb);
        obj.ub0 = fastopt.flattenstruct(arg.ub);

        % Initialize state.
        obj.values = obj.init0;
        obj.lb = obj.lb0;
        obj.ub = obj.ub0;
        obj.startFlag = [];

        % Initialize the GUI.
        [obj.freeParam, obj.dof, obj.mult] = obj.getFreeParams(arg.modelspec);
        obj.TdegC = arg.modelspec.temps-273.15;
        obj.buildToplevelGUI();
        obj.buildControlPanel();
        obj.buildParameterPanel();
        obj.initializePlot();
        obj.updatePlot();
    end % constructor

    function data = run(obj)
        % Run main GUI loop.
        obj.status.Text = "Idle";
        opttype = [];
        while ishandle(obj.fig)
            pause(0.01);
            drawnow;
            if ~isempty(obj.startFlag)
                opttype = obj.startFlag;
                obj.startFlag = [];
                break;
            end
        end % while

        data = struct;
        data.init = fastopt.unflattenstruct(obj.values);
        data.lb = fastopt.unflattenstruct(obj.lb);
        data.ub = fastopt.unflattenstruct(obj.ub);
        if ~ishandle(obj.fig)
            % User canceled the optimization.
            data.OptimizationType = [];
            data.SwarmSize = [];
            data.SwarmMaximumIterations = [];
            data.FminMaximumIterations = [];
            data.FminMaximumFunctionEvals = [];
            data.updateControlsFcn = [];
            data.updateStatusFcn = [];
            data.stop = true;
        else
            data.OptimizationType = opttype;
            data.SwarmSize = round(abs(obj.edt.popsize.Value));
            data.SwarmMaximumIterations = round(abs(obj.edt.maxiter.Value));
            data.FminMaximumIterations = round(abs(obj.edt.fminmaxiter.Value));
            data.FminMaximumFunctionEvals = round(abs(obj.edt.fminmaxevals.Value));
            [data.updateControlsFcn, data.updateStatusFcn] = obj.updateGUIFactory(obj);
            data.stop = false;
        end
    end % run()

    function update(obj, values)
        obj.values = fastopt.flattenstruct(values);
        obj.updateParameterGUIFromState();
        obj.updatePlot();
    end

    function hide(obj)
        obj.fig.Visible = 'off';
        drawnow;
    end

    function show(obj)
        obj.fig.Visible = 'on';
        drawnow;
    end

    function lock(obj, stopFcn)
        edtnames = fieldnames(obj.edt);
        for k = 1:length(edtnames)
            edtname = edtnames{k};
            edtfield = obj.edt.(edtname);
            edtfield.Enable = 'off';
            edtfield.Editable = 'off';
        end
        sldnames = fieldnames(obj.sld);
        for k = 1:length(sldnames)
            sldname = sldnames{k};
            sldfield = obj.sld.(sldname);
            sldfield.Enable = 'off';
        end
        btnnames = fieldnames(obj.btn);
        for k = 1:length(btnnames)
            btnname = btnnames{k};
            btnfield = obj.btn.(btnname);
            btnfield.Enable = 'off';
        end
        % Set stop button callback.
        obj.btn.stop.ButtonPushedFcn = @stopFcnWrapper;
        obj.btn.stop.Enable = 'on';
        drawnow;
        function stopFcnWrapper(src,~)
            src.Enable = 'off';
            drawnow;
            stopFcn();
        end
    end

    function unlock(obj)
        edtnames = fieldnames(obj.edt);
        for k = 1:length(edtnames)
            edtname = edtnames{k};
            edtfield = obj.edt.(edtname);
            edtfield.Enable = 'on';
            edtfield.Editable = 'on';
        end
        sldnames = fieldnames(obj.sld);
        for k = 1:length(sldnames)
            sldname = sldnames{k};
            sldfield = obj.sld.(sldname);
            sldfield.Enable = 'on';
        end
        btnnames = fieldnames(obj.btn);
        for k = 1:length(btnnames)
            btnname = btnnames{k};
            btnfield = obj.btn.(btnname);
            btnfield.Enable = 'on';
        end
        % Clear stop button callback.
        obj.btn.ButtonPushedFcn = []; 
        obj.btn.stop.Enable = 'off';
        drawnow;
    end
end

methods(Access=protected)
    function buildToplevelGUI(obj)
        obj.fig = uifigure("Name",obj.arg.WindowName);
        obj.fig.Units = 'normalized';
        obj.fig.Position = [0.05 0.05 0.9 0.9];
        grid = uigridlayout(obj.fig,[2 2]);
        grid.ColumnWidth = {'1x','2x'};
        grid.RowHeight = {150,'1x'};
        obj.panel.control = uipanel(grid);
        obj.panel.control.Layout.Row = 1;
        obj.panel.control.Layout.Column = 1;
        obj.panel.parameter = uipanel(grid);
        obj.panel.parameter.Layout.Row = 2;
        obj.panel.parameter.Layout.Column = 1;
        obj.panel.plot = uipanel(grid);
        obj.panel.plot.Layout.Row = [1 2];
        obj.panel.plot.Layout.Column = 2;
    end % buildToplevelGUI()

    function buildControlPanel(obj)
        grid = uigridlayout(obj.panel.control,[4 6]);
        grid.RowHeight = {'1.5x','1.5x','1.5x','1x'};
        obj.btn.load = uibutton(grid, ...
            "Text","Load Parameter Values", ...
            "ButtonPushedFcn",@obj.loadButtonPushed);
        obj.btn.load.Layout.Column = [1 2];
        obj.btn.save = uibutton(grid, ...
            "Text","Save Parameter Values", ...
            "ButtonPushedFcn",@obj.saveButtonPushed);
        obj.btn.save.Layout.Column = [3 4];
        obj.btn.reset = uibutton(grid, ...
            "Text","Reset Parameter Values", ...
            "ButtonPushedFcn",@obj.resetButtonPushed);
        obj.btn.reset.Layout.Column = [5 6];
        obj.lab.popsize = uilabel(grid,"Text","Swarm Size");
        obj.edt.popsize = uieditfield(grid,'numeric', ...
            "Value",obj.arg.SwarmSize);
        obj.lab.maxiter = uilabel(grid,"Text","Max. Iterations");
        obj.edt.maxiter = uieditfield(grid,'numeric', ...
            "Value",obj.arg.SwarmMaximumIterations);
        obj.btn.start = uibutton(grid, ...
            "Text","Start", ...
            "ButtonPushedFcn",@obj.startButtonPushed);
        obj.btn.stop = uibutton(grid, ...
            "Text","Stop","Enable","off");
        obj.btn.stop.Layout.Row = [2 3];
        obj.btn.stop.Layout.Column = 6;
        obj.lab.fminmaxiter = uilabel(grid,"Text","Fmin Max. Iter.");
        obj.lab.fminmaxiter.Layout.Row = 3;
        obj.lab.fminmaxiter.Layout.Column = 1;
        obj.edt.fminmaxiter = uieditfield(grid,'numeric', ...
            "Value",obj.arg.FminMaximumIterations);
        obj.edt.fminmaxiter.Layout.Row = 3;
        obj.edt.fminmaxiter.Layout.Column = 2;
        obj.lab.fminmaxevals = uilabel(grid,"Text","Max. F-Evals");
        obj.lab.fminmaxevals.Layout.Row = 3;
        obj.lab.fminmaxevals.Layout.Column = 3;
        obj.edt.fminmaxevals = uieditfield(grid,'numeric', ...
            "Value",obj.arg.FminMaximumFunctionEvals);
        obj.edt.fminmaxevals.Layout.Row = 3;
        obj.edt.fminmaxevals.Layout.Column = 4;
        obj.btn.start2 = uibutton(grid, ...
            "Text","Start", ...
            "ButtonPushedFcn",@obj.start2ButtonPushed);
        obj.btn.start2.Layout.Row = 3;
        obj.btn.start2.Layout.Column = 5;
        obj.status = uilabel(grid,"Text","Idle");
        obj.status.Layout.Column = [1 6];
    end % buildControlPanel()

    function buildParameterPanel(obj)
        grid = uigridlayout(obj.panel.parameter,[obj.dof+1 5]);
        h1 = uilabel(grid,"Text","Param.","FontWeight","bold");
        h2 = uilabel(grid,"Text","Value","FontWeight","bold");
        h2.Layout.Column = [2 3];
        h3 = uilabel(grid,"Text","Lower","FontWeight","bold");
        h4 = uilabel(grid,"Text","Upper","FontWeight","bold");

        % Iterate free parameters.
        cursor = 1;
        paramNames = fieldnames(obj.freeParam);
        for indParam = 1:length(paramNames)
            paramName = paramNames{indParam};
            paramNameDisplay = strrep(paramName,'__','.');
            param = obj.freeParam.(paramName);

            % Determine multiplicity.
            multiplicity = 1;
            if strcmpi(param.tempfcn,'lut')
                multiplicity = obj.mult;  % need fields for each temperature
            end

            % Check for activation energy, set-up loop iterations.
            iterations = multiplicity;
            if strcmpi(param.tempfcn,'Eact')
                % Add extra iteration for activation-energy field.
                iterations = iterations + 1;
            end
        
            % Iterate degrees of freedom within free parameter.
            for indItem = 1:param.len
                % Iterate temperature multiplicity of each item.
                for m = 1:iterations
                    if m <= multiplicity
                        pname = paramName;
                        % Iterations for paramter at each temperature.
                        itemName = sprintf('%s_%d_%d',paramName,indItem,m);
                        if param.len > 1
                            itemNameDisplay = sprintf('%s%d', ...
                                paramNameDisplay,indItem);
                        else
                            itemNameDisplay = paramNameDisplay;
                        end
                        if multiplicity>1
                            itemNameDisplay = sprintf('%s @%.1fC', ...
                                itemNameDisplay,obj.TdegC(m));
                        end
                        if param.logscale
                            itemNameDisplay = [itemNameDisplay ' (log)'];
                        end
                        data = struct;
                        data.nameDisplay = itemNameDisplay;
                        data.param = param;
                        data.paramName = paramName;
                        data.indParam = indParam;
                        data.indItem = indItem;
                        data.indMult = m;
                    else
                        if indItem ~= param.len
                            continue;
                        end
                        % Extra iteration for activation energy field.
                        % Always use linear scale.
                        pname = [paramName '_Eact'];
                        itemName = [paramName '_Eact'];
                        itemNameDisplay = [paramNameDisplay '_Eact'];
                        data = struct;
                        data.nameDisplay = itemNameDisplay;
                        data.param = param;
                        data.param.logscale = false;
                        data.paramName = itemName;
                        data.indParam = indParam;
                        data.indItem = 1;
                        data.indMult = 1;
                    end
                    data.labName = itemName;
                    data.sldName = itemName;
                    data.edtValueName = [itemName '__value'];
                    data.edtLowerName = [itemName '__lower'];
                    data.edtUpperName = [itemName '__upper'];

                    % Fetch values for input fields.
                    initial = obj.init0.(pname);
                    lower = obj.lb0.(pname);
                    upper = obj.ub0.(pname);
                    initial(initial<lower) = lower(initial<lower);
                    initial(initial>upper) = upper(initial>upper);
                    initialPct = obj.toSliderValue( ...
                        initial,lower,upper,data.param.logscale);
            
                    % Parameter label.
                    obj.lab.(data.labName) = uilabel( ...
                        grid,"Text",itemNameDisplay);
                    obj.lab.(data.labName).Layout.Row = cursor+1;
                    obj.lab.(data.labName).Layout.Column = 1;
                
                    % Slider for adjusting initial value.
                    obj.sld.(data.sldName) = uislider( ...
                        grid,"Value",initialPct(data.indItem));
                    obj.sld.(data.sldName).Layout.Row = cursor+1;
                    obj.sld.(data.sldName).Layout.Column = 2;
                    obj.sld.(data.sldName).UserData = data;
                    obj.sld.(data.sldName).UserData.type = "value";
                    obj.sld.(data.sldName).ValueChangedFcn = @obj.updateSlider;
                    obj.sld.(data.sldName).MajorTickLabels = {};
            
                    % Initial value editor.
                    obj.edt.(data.edtValueName) = uieditfield( ...
                        grid,'numeric',"Value",initial(data.indItem));
                    obj.edt.(data.edtValueName).Layout.Row = cursor+1;
                    obj.edt.(data.edtValueName).Layout.Column = 3;
                    obj.edt.(data.edtValueName).UserData = data;
                    obj.edt.(data.edtValueName).UserData.type = "value";
                    obj.edt.(data.edtValueName).ValueChangedFcn = @obj.updateEditValue;
            
                    % Lower bound editor.
                    obj.edt.(data.edtLowerName) = uieditfield( ...
                        grid,'numeric',"Value",lower(data.indItem));
                    obj.edt.(data.edtLowerName).Layout.Row = cursor+1;
                    obj.edt.(data.edtLowerName).Layout.Column = 4;
                    obj.edt.(data.edtLowerName).UserData = data;
                    obj.edt.(data.edtLowerName).UserData.type = "lower";
                    obj.edt.(data.edtLowerName).ValueChangedFcn = @obj.updateEditLower;
            
                    % Upper bound editor.
                    obj.edt.(data.edtUpperName) = uieditfield( ...
                        grid,'numeric',"Value",upper(data.indItem));
                    obj.edt.(data.edtUpperName).Layout.Row = cursor+1;
                    obj.edt.(data.edtUpperName).Layout.Column = 5;
                    obj.edt.(data.edtUpperName).UserData = data;
                    obj.edt.(data.edtUpperName).UserData.type = "upper";
                    obj.edt.(data.edtUpperName).ValueChangedFcn = @obj.updateEditUpper;
            
                    cursor = cursor + 1;
                end % for mult
            end % for item
        end % for param

        % Do the scrolling after component creation, or it will take 
        % forever!
        grid.ColumnWidth = {'1.5x','2x','1x','1x','1x'};
        grid.Scrollable = 'on';
        grid.RowHeight = num2cell(20*ones(1,cursor));
    end % buildParameterPanel()

    function initializePlot(obj)
        obj.arg.PlotInitializeFcn(obj.panel.plot);
    end
    
    function updatePlot(obj)
        obj.arg.PlotUpdateFcn(fastopt.unflattenstruct(obj.values));
    end

    function updateParameterGUIFromState(obj)
        % Update value, lower, upper fields.
        edtnames = fieldnames(obj.edt);
        for k = 1:length(edtnames)
            edtname = edtnames{k};
            edtfield = obj.edt.(edtname);
            if ~isstruct(edtfield.UserData) || ~isfield(edtfield.UserData,'type')
                % Non-parameter edit field.
                continue;
            end
            data = edtfield.UserData;
            if strcmpi(data.type,"value")
                value = obj.values.(data.paramName)(data.indItem,data.indMult);
                obj.edt.(edtname).Value = value;
            elseif strcmpi(data.type,"lower")
                lwr = obj.lb.(data.paramName)(data.indItem,data.indMult);
                obj.edt.(edtname).Value = lwr;
            elseif strcmpi(data.type,"upper")
                upr = obj.ub.(data.paramName)(data.indItem,data.indMult);
                obj.edt.(edtname).Value = upr;
            end
        end % for

        % Update sliders.
        sldnames = fieldnames(obj.sld);
        for k = 1:length(sldnames)
            sldname = sldnames{k};
            sldfield = obj.sld.(sldname);
            if ~isstruct(sldfield.UserData)
                % Non-parameter slider field.
                continue;
            end
            data = sldfield.UserData;
            lwr = obj.lb.(data.paramName)(data.indItem,data.indMult);
            upr = obj.ub.(data.paramName)(data.indItem,data.indMult);
            value = obj.values.(data.paramName)(data.indItem,data.indMult);
            valuePct = obj.toSliderValue(value,lwr,upr,data.param.logscale);
            obj.sld.(sldname).Value = valuePct;
        end % for
    end % updateParameterValues()

    % Event handlers ------------------------------------------------------
    function startButtonPushed(obj,~,~)
        obj.startFlag = 'particleswarm';
    end

    function start2ButtonPushed(obj,~,~)
        obj.startFlag = 'fmincon';
    end

    function stopButtonPushed(obj,~,~)
        obj.stopFlag = true;
    end

    function loadButtonPushed(obj,~,~)
        [filename,filepath] = uigetfile("*.mat");
    
        % Prevents ui window from falling into background after uigetfile.
        drawnow;
        figure(obj.fig);
    
        if isnumeric(filename) && filename == 0
            % User canceled the operation.
            return;
        end
        dat = load([filepath,filename]);
        obj.lb = fastopt.flattenstruct(dat.lb);
        obj.ub = fastopt.flattenstruct(dat.ub);
        obj.values = fastopt.flattenstruct(dat.values);
        obj.updateParameterGUIFromState();
        obj.updatePlot();
    end
    
    function saveButtonPushed(obj,~,~)
        [filename,filepath] = uiputfile("*.mat");
    
        % Prevents ui window from falling into background after uiputfile.
        drawnow;
        figure(obj.fig);
    
        if isnumeric(filename) && filename == 0
            % User canceled the operation.
            return;
        end
        dat.lb = fastopt.unflattenstruct(obj.lb);
        dat.ub = fastopt.unflattenstruct(obj.ub);
        dat.values = fastopt.unflattenstruct(obj.values);
        save([filepath,filename],'-struct','dat');
    end
    
    function resetButtonPushed(obj,~,~)
        obj.lb = obj.lb0;
        obj.ub = obj.ub0;
        obj.values = obj.init0;
        obj.updateParameterGUIFromState();
        obj.updatePlot();
    end
    
    function updateSlider(obj,src,event)
        data = src.UserData;
    
        % Compute new value from slider setting.
        lwr = obj.lb.(data.paramName)(data.indItem,data.indMult);
        upr = obj.ub.(data.paramName)(data.indItem,data.indMult);
        newValuePct = event.Value;
        newValue = obj.fromSliderValue(newValuePct,lwr,upr,data.param.logscale);
    
        % Update state.
        obj.values.(data.paramName)(data.indItem,data.indMult) = newValue;
    
        % Update GUI components.
        obj.edt.(data.edtValueName).Value = newValue;
        obj.updatePlot();
    end
    
    function updateEditValue(obj,src,event)
        data = src.UserData;
    
        % Compute slider setting from new initial value.
        lwr = obj.lb.(data.paramName)(data.indItem,data.indMult);
        upr = obj.ub.(data.paramName)(data.indItem,data.indMult);
        newValue = event.Value;
        if newValue < lwr
            newValue = lwr;   % ensure greater than lower bound
        elseif newValue > upr
            newValue = upr;   % ensure less than upper bound
        end
        newValuePct = obj.toSliderValue(newValue,lwr,upr,data.param.logscale);
    
        % Update state.
        obj.values.(data.paramName)(data.indItem,data.indMult) = newValue;
    
        % Update GUI components.
        obj.edt.(data.edtValueName).Value = newValue;  % in case value was clamped
        obj.sld.(data.sldName).Value = newValuePct;
        obj.updatePlot();
    end
    
    function updateEditLower(obj,src,event)
        data = src.UserData;
    
        % Compute slider setting from new lower bound.
        newLwr = event.Value;
        upr = obj.ub.(data.paramName)(data.indItem,data.indMult);
        value = obj.values.(data.paramName)(data.indItem,data.indMult);
        if value < newLwr
            value = newLwr;  % ensure greater than lower bound
        end
        valuePct = obj.toSliderValue(value,newLwr,upr,data.param.logscale);
    
        % Update state.
        obj.lb.(data.paramName)(data.indItem,data.indMult) = newLwr;
        obj.values.(data.paramName)(data.indItem,data.indMult) = value; % in case value was clamped
    
        % Update GUI components.
        obj.edt.(data.edtValueName).Value = value;
        obj.sld.(data.sldName).Value = valuePct;
        obj.updatePlot();
    end
    
    function updateEditUpper(obj,src,event)
        data = src.UserData;
    
        % Compute slider setting from new lower bound.
        newUpr = event.Value;
        lwr = obj.lb.(data.paramName)(data.indItem,data.indMult);
        value = obj.values.(data.paramName)(data.indItem,data.indMult);
        if value > newUpr
            value = newUpr;  % ensure less than lower bound
        end
        valuePct = obj.toSliderValue(value,lwr,newUpr,data.param.logscale);
    
        % Update upper bound.
        obj.ub.(data.paramName)(data.indItem,data.indMult) = newUpr;
        obj.values.(data.paramName)(data.indItem,data.indMult) = value;  % in case value was clamped
    
        % Update GUI components.
        obj.edt.(data.edtValueName).Value = value;
        obj.sld.(data.sldName).Value = valuePct;
        obj.updatePlot();
    end
end

methods(Static,Access=protected)
    function [updateControlsFcn, updateStatusFcn] = updateGUIFactory(obj)
        edt = obj.edt;
        sld = obj.sld;
        status = obj.status;
        lb = obj.lb;
        ub = obj.ub;
        toSliderValue = @obj.toSliderValue;
        updateControlsFcn = @updateControls;
        updateStatusFcn = @updateStatus;

        function updateStatus(statusString)
            status.Text = statusString;
            drawnow;
        end

        function updateControls(newValues)
            newValues = fastopt.flattenstruct(newValues);

            % Update value fields.
            edtnames = fieldnames(edt);
            for k = 1:length(edtnames)
                edtname = edtnames{k};
                edtfield = edt.(edtname);
                if ~isstruct(edtfield.UserData) || ...
                   ~isfield(edtfield.UserData,'type')
                    % Non-parameter edit field.
                    continue;
                end
                data = edtfield.UserData;
                if strcmpi(data.type,"value")
                    value = newValues.(data.paramName)(data.indItem,data.indMult);
                    edt.(edtname).Value = value;
                end
            end % for
    
            % Update sliders.
            sldnames = fieldnames(sld);
            for k = 1:length(sldnames)
                sldname = sldnames{k};
                sldfield = obj.sld.(sldname);
                if ~isstruct(sldfield.UserData)
                    % Non-parameter slider field.
                    continue;
                end
                data = sldfield.UserData;
                lwr = lb.(data.paramName)(data.indItem,data.indMult);
                upr = ub.(data.paramName)(data.indItem,data.indMult);
                value = newValues.(data.paramName)(data.indItem,data.indMult);
                valuePct = toSliderValue(value,lwr,upr,data.param.logscale);
                sld.(sldname).Value = valuePct;
            end % for

            drawnow;
        end % updateParameterValues()
    end

    function [params, dof, mult] = getFreeParams(modelspec)
        %GETFREEPARAMS Get parameters that are allowed to vary in the
        %  optimization and the number of degrees of freedom (may differ
        %  from number of free params b/c some params may have length > 1).
        params = struct;
        paramNames = fieldnames(modelspec.params);
        for indItem = 1:length(paramNames)
            paramName = paramNames{indItem};
            param = modelspec.params.(paramName);
            if ~isfield(param,'fix')
                params.(paramName) = param;
            end
        end % for
        dof = modelspec.nvars;
        mult = modelspec.ntemps;
    end % getFreeParams()

    function value = fromSliderValue(valueSlider,lwr,upr,logscale)
        %FROMSLIDERVALUE Convert percent value from slider to absolute units.
        if logscale
            value = log10(lwr)+(valueSlider/100).*(log10(upr)-log10(lwr));
            value = 10^value;
        else
            value = lwr+(valueSlider/100).*(upr-lwr);
        end
    end % fromSliderValue()
    
    function valueSlider = toSliderValue(value,lwr,upr,logscale)
        %TOSLIDERVALUE Convert absolute value to percent value for slider.
        ind = lwr==upr;
        if any(ind)
            upr(ind) = lwr(ind)+1;
        end
        if logscale
            valueSlider = ... 
                100*(log10(value)-log10(lwr))./(log10(upr)-log10(lwr));
        else
            valueSlider = 100*(value-lwr)./(upr-lwr);
        end
    end % toSliderValue()
end
end % classdef