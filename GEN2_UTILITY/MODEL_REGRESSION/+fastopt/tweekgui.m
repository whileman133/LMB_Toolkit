function [init,lb,ub] = tweekgui(modelspec,init,lb,ub,initializefcn,updatefcn)
%TWEEKGUI Construct a GUI for adjusting initial parameter values / bounds.

init = fastopt.flattenstruct(init);
lb = fastopt.flattenstruct(lb);
ub = fastopt.flattenstruct(ub);

% First, fetch all free parameters.
params = struct;
paramNames = fieldnames(modelspec.params);
nparam = 0;
for indItem = 1:length(paramNames)
    paramName = paramNames{indItem};
    param = modelspec.params.(paramName);
    if ~isfield(param,'fix')
        params.(paramName) = param;
        nparam = nparam + param.len;
    end
end % for
paramNames = fieldnames(params);

% Construct figure window.
fig = uifigure("Name","Optimization Tweeker");
fig.Units = 'normalized';
fig.Position = [0.05 0.05 0.9 0.9];
gridtop = uigridlayout(fig,[1 2]);
gridtop.ColumnWidth = {'1x','2x'};
panelparam = uipanel(gridtop);
paneldisplay = uipanel(gridtop);

% Construct parameter panel for adjusting init,lb,ub.
gridparam = uigridlayout(panelparam,[nparam+1 5]);
gridparam.ColumnWidth = {'1x','2x','1x','1x','1x'};
h1 = uilabel(gridparam,"Text","Param.","FontWeight","bold");
h2 = uilabel(gridparam,"Text","Value","FontWeight","bold");
h2.Layout.Column = [2 3];
h3 = uilabel(gridparam,"Text","Lower","FontWeight","bold");
h4 = uilabel(gridparam,"Text","Upper","FontWeight","bold");
cursor = 1;
sldinit = gobjects(nparam,1);
editinit = gobjects(nparam,1);
editlower = gobjects(nparam,1);
editupper = gobjects(nparam,1);
for indParam = 1:length(paramNames)
    paramName = paramNames{indParam};
    paramNameDisplay = strrep(paramName,'__','.');
    param = params.(paramName);
    initial = init.(paramName);
    lower = lb.(paramName);
    upper = ub.(paramName);
    initialPct = toSliderValue(initial,lower,upper,param.logscale);

    for indItem = 1:param.len
        if param.len > 1
            itemNameDisplay = sprintf('%s%d',paramNameDisplay,indItem);
        else
            itemNameDisplay = paramNameDisplay;
        end
        data = struct;
        data.nameDisplay = itemNameDisplay;
        data.param = param;
        data.paramName = paramName;
        data.indParam = indParam;
        data.indItem = indItem;
        data.indGraphics = cursor;

        % Parameter label.
        lab = uilabel( ...
            gridparam,"Text",itemNameDisplay);
        lab.Layout.Row = cursor+1;
        lab.Layout.Column = 1;
    
        % Slider for adjusting initial value.
        sldinit(cursor) = uislider( ...
            gridparam,"Value",initialPct(indItem));
        sldinit(cursor).Layout.Row = cursor+1;
        sldinit(cursor).Layout.Column = 2;
        sldinit(cursor).UserData = data;
        sldinit(cursor).ValueChangedFcn = @updateSlider;

        % Initial value editor.
        editinit(cursor) = uieditfield( ...
            gridparam,'numeric',"Value",initial(indItem));
        editinit(cursor).Layout.Row = cursor+1;
        editinit(cursor).Layout.Column = 3;
        editinit(cursor).UserData = data;
        editinit(cursor).ValueChangedFcn = @updateEditInitial;

        % Lower bound editor.
        editlower(cursor) = uieditfield( ...
            gridparam,'numeric',"Value",lower(indItem));
        editlower(cursor).Layout.Row = cursor+1;
        editlower(cursor).Layout.Column = 4;
        editlower(cursor).UserData = data;
        editlower(cursor).ValueChangedFcn = @updateEditLower;

        % Upper bound editor.
        editupper(cursor) = uieditfield( ...
            gridparam,'numeric',"Value",upper(indItem));
        editupper(cursor).Layout.Row = cursor+1;
        editupper(cursor).Layout.Column = 5;
        editupper(cursor).UserData = data;
        editupper(cursor).ValueChangedFcn = @updateEditUpper;

        cursor = cursor + 1;
    end % for
end

% Construct display panel.
initializeDisplay();
updateDisplay();

% Wait until user finishes adjusting parameter values.
drawnow;
uiwait(fig);

% Unflatten init,lb,ub structures.
init = fastopt.unflattenstruct(init);
lb = fastopt.unflattenstruct(lb);
ub = fastopt.unflattenstruct(ub);

% Event handlers ----------------------------------------------------------
function updateSlider(src,event)
    p = src.UserData.param;
    pname = src.UserData.paramName;
    indI = src.UserData.indItem;
    indG = src.UserData.indGraphics;

    % Compute new value from slider setting.
    lwr = lb.(pname)(indI);
    upr = ub.(pname)(indI);
    newValuePct = event.Value;
    newInit = fromSliderValue(newValuePct,lwr,upr,p.logscale);

    % Update initial value.
    init.(pname)(indI) = newInit;

    % Update GUI components.
    editinit(indG).Value = newInit;
    updateDisplay();
end

function updateEditInitial(src,event)
    p = src.UserData.param;
    pname = src.UserData.paramName;
    indI = src.UserData.indItem;
    indG = src.UserData.indGraphics;

    % Compute slider setting from new initial value.
    lwr = lb.(pname)(indI);
    upr = ub.(pname)(indI);
    newInit = event.Value;
    if newInit < lwr
        newInit = lwr;   % ensure greater than lower bound
    elseif newInit > upr
        newInit = upr;   % ensure less than upper bound
    end
    newInitPct = toSliderValue(newInit,lwr,upr,p.logscale);

    % Update initial value.
    init.(pname)(indI) = newInit;
    editinit(indG).Value = newInit;   % in case value was clamped

    % Update GUI components.
    sldinit(indG).Value = newInitPct;
    updateDisplay();
end

function updateEditLower(src,event)
    p = src.UserData.param;
    pname = src.UserData.paramName;
    indI = src.UserData.indItem;
    indG = src.UserData.indGraphics;

    % Compute slider setting from new lower bound.
    newLwr = event.Value;
    upr = ub.(pname)(indI);
    newInit = init.(pname)(indI);
    if newInit < newLwr
        newInit = newLwr;  % ensure greater than lower bound
    end
    newInitPct = toSliderValue(newInit,newLwr,upr,p.logscale);

    % Update lower bound.
    lb.(pname)(indI) = newLwr;
    init.(pname)(indI) = newInit;    % in case value was clamped

    % Update GUI components.
    editinit(indG).Value = newInit;
    sldinit(indG).Value = newInitPct;
    updateDisplay();
end

function updateEditUpper(src,event)
    p = src.UserData.param;
    pname = src.UserData.paramName;
    indI = src.UserData.indItem;
    indG = src.UserData.indGraphics;

    % Compute slider setting from new lower bound.
    newUpr = event.Value;
    lwr = lb.(pname)(indI);
    newInit = init.(pname)(indI);
    if newInit > newUpr
        newInit = newUpr;  % ensure less than lower bound
    end
    newInitPct = toSliderValue(newInit,lwr,newUpr,p.logscale);

    % Update upper bound.
    ub.(pname)(indI) = newUpr;
    init.(pname)(indI) = newInit;    % in case value was clamped

    % Update GUI components.
    editinit(indG).Value = newInit;
    sldinit(indG).Value = newInitPct;
    updateDisplay();
end

% Utility functions -------------------------------------------------------
function initializeDisplay()
    initializefcn( ...
        paneldisplay, ...
        fastopt.unflattenstruct(init), ...
        fastopt.unflattenstruct(lb), ...
        fastopt.unflattenstruct(ub));
end

function updateDisplay()
    updatefcn( ...
        paneldisplay, ...
        fastopt.unflattenstruct(init), ...
        fastopt.unflattenstruct(lb), ...
        fastopt.unflattenstruct(ub));
end

function value = fromSliderValue(valueSlider,lwr,upr,logscale)
    %FROMSLIDERVALUE Convert percent value from slider to absolute units.
    if logscale
        value = log10(lwr)+(valueSlider/100).*(log10(upr)-log10(lwr));
        value = 10^value;
    else
        value = lwr+(valueSlider/100).*(upr-lwr);
    end
end

function valueSlider = toSliderValue(value,lwr,upr,logscale)
    %TOSLIDERVALUE Convert absolute value to percent value for slider.
    if logscale
        valueSlider = ... 
            100*(log10(value)-log10(lwr))./(log10(upr)-log10(lwr));
    else
        valueSlider = 100*(value-lwr)./(upr-lwr);
    end
end

end