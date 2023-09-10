function [initializeFcn, updateFcn] = uiH2ImpedanceCallbacks(modelspec,Zref,freqRef,socPctRef,TdegC)
%UIH2IMPEDANCECALLBACKS

nsoc = length(socPctRef{1});
J = modelspec.params.pos__U0.len;

% Define globals.
axZ = gobjects(nsoc,1);
zplotselect = gobjects(1,1);
tempselect = gobjects(1,1);
gridz = gobjects(1,1);
linesZ = gobjects(nsoc,1);
lastParamValues = [];
[~,tempSelect] = min(abs(TdegC-25));  % choose temp closest to 25degC by default.
plotSelect = "Nyq";

initializeFcn = @initializeUIFig;
updateFcn = @updateUIFig;

function initializeUIFig(parent)
    % Construct layout.
    gridtop = uigridlayout(parent,[2 2]);
    gridtop.ColumnWidth = {'3.5x','1x'};
    gridtop.RowHeight = {'1x', '3x'};
    panelz = uipanel(gridtop,'BorderType','none');
    panelz.Layout.Row = [1 2];
    panelz.Layout.Column = 1;
    panelcontrols = uipanel(gridtop,'BorderType','none');
    panelcontrols.Layout.Row = 1;
    panelcontrols.Layout.Column = 2;
    panelct = uipanel(gridtop,'BorderType','none');
    panelct.Layout.Row = 2;
    panelct.Layout.Column = 2;

    % Construct controls.
    gridcontrols = uigridlayout(panelcontrols,[3 1]);
    gridcontrols.RowHeight = {30, 30, '1x'};
    zplotselect = uidropdown(gridcontrols, ...
        "Items",["Nyquist","Bode Magnitude","Bode Phase","Real Part","Imaginary Part"], ...
        "ItemsData",["Nyq","BodeMag","BodePhase","Real","Imag"],...
        "ValueChangedFcn",@(src,event)updateZPlotType(event.Value));
    tempselect = uidropdown(gridcontrols, ...
        "Items",arrayfun(@(x)sprintf("T=%.1fdegC",x),TdegC), ...
        "ItemsData",1:length(TdegC),"Value",tempSelect,...
        "ValueChangedFcn",@(src,event)updateTempSelect(event.Value));

    % Construct impedance plots.
    gridz = uigridlayout(panelz,[5 ceil(nsoc/5)]);
    for k = 1:nsoc
        ax = uiaxes(gridz);
        axZ(k) = ax;
        formatAxes(ax);
        plot(ax,real(Zref{tempSelect}(:,k)),-imag(Zref{tempSelect}(:,k)),'b.');
        hold(ax,'on');
        linesZ(k) = plot(ax,NaN,NaN,'r-');
        title(ax,sprintf('%.0f%% SOC',socPctRef{tempSelect}(k)));
        setAxesNyquist('axes',ax);
    end
end % initializeUIPlot()

function updateUIFig(paramValues)
    % Save parameter values so we can update the plots when plot type
    % changes.
    lastParamValues = paramValues;

    % Select appropriate data-set.
    paramValues = fastopt.splittemps(paramValues,modelspec);
    paramValues = paramValues(tempSelect);

    % Calculate impedance predicted by the linear EIS model.
    Zmodel = getH2Impedance( ...
        paramValues,freqRef{tempSelect},socPctRef{tempSelect},TdegC(tempSelect));

    % Update model predictions on plots.
    if plotSelect == "Nyq"
        for idxSOC = 1:nsoc
            linesZ(idxSOC).XData = real(Zmodel(:,idxSOC));
            linesZ(idxSOC).YData = -imag(Zmodel(:,idxSOC));
        end
    elseif plotSelect == "BodeMag"
        for idxSOC = 1:nsoc
            linesZ(idxSOC).XData = freqRef{tempSelect};
            linesZ(idxSOC).YData = abs(Zmodel(:,idxSOC));
        end
    elseif plotSelect == "BodePhase"
        for idxSOC = 1:nsoc
            linesZ(idxSOC).XData = freqRef{tempSelect};
            linesZ(idxSOC).YData = angle(Zmodel(:,idxSOC))*180/pi;
        end
    elseif plotSelect == "Real"
        for idxSOC = 1:nsoc
            linesZ(idxSOC).XData = freqRef{tempSelect};
            linesZ(idxSOC).YData = real(Zmodel(:,idxSOC));
        end
    elseif plotSelect == "Imag"
        for idxSOC = 1:nsoc
            linesZ(idxSOC).XData = freqRef{tempSelect};
            linesZ(idxSOC).YData = imag(Zmodel(:,idxSOC));
        end
    end
end % updateUIPlot()

function updateZPlotType(plottype)
    plotSelect = plottype;
    redrawPlots();
end % updateZPlotType()

function updateTempSelect(indT)
    tempSelect = indT;
    redrawPlots();
end % updateZPlotType()

function redrawPlots()
    if plotSelect == "Nyq"    
        for k = 1:nsoc
            ax = axZ(k);
            cla(ax);
            ax.XScale = 'linear';
            ax.YScale = 'linear';
            plot(ax,real(Zref{tempSelect}(:,k)),-imag(Zref{tempSelect}(:,k)),'b.');
            hold(ax,'on');
            linesZ(k) = plot(ax,NaN,NaN,'r-');
            title(ax,sprintf('%.0f%% SOC',socPctRef{tempSelect}(k)));
            setAxesNyquist('axes',ax);
        end 
    elseif plotSelect == "BodeMag"
        for k = 1:nsoc
            ax = axZ(k);
            cla(ax);
            ax.XScale = 'log';
            ax.YScale = 'log';
            loglog(ax,freqRef{tempSelect},abs(Zref{tempSelect}(:,k)),'b.');
            hold(ax,'on');
            linesZ(k) = loglog(ax,NaN,NaN,'r-');
            title(ax,sprintf('%.0f%% SOC',socPctRef{tempSelect}(k)));
            ylim(ax,[min(abs(Zref{tempSelect}(:,k))) max(abs(Zref{tempSelect}(:,k)))]);
            xlim(ax,[min(freqRef{tempSelect}) max(freqRef{tempSelect})]);
        end 
    elseif plotSelect == "BodePhase"
        for k = 1:nsoc
            ax = axZ(k);
            cla(ax);
            ax.XScale = 'log';
            ax.YScale = 'linear';
            semilogx(ax,freqRef{tempSelect},angle(Zref{tempSelect}(:,k))*180/pi,'b.');
            hold(ax,'on');
            linesZ(k) = semilogx(ax,NaN,NaN,'r-');
            title(ax,sprintf('%.0f%% SOC',socPctRef{tempSelect}(k)));
            ylim(ax,[min(angle(Zref{tempSelect}(:,k))*180/pi) max(angle(Zref{tempSelect}(:,k))*180/pi)]);
            xlim(ax,[min(freqRef{tempSelect}) max(freqRef{tempSelect})]);
        end
    elseif plotSelect == "Real"
        for k = 1:nsoc
            ax = axZ(k);
            cla(ax);
            ax.XScale = 'log';
            ax.YScale = 'linear';
            loglog(ax,freqRef{tempSelect},real(Zref{tempSelect}(:,k)),'b.');
            hold(ax,'on');
            linesZ(k) = semilogx(ax,NaN,NaN,'r-');
            title(ax,sprintf('%.0f%% SOC',socPctRef{tempSelect}(k)));
            ylim(ax,[min(real(Zref{tempSelect}(:,k))) max(real(Zref{tempSelect}(:,k)))]);
            xlim(ax,[min(freqRef{tempSelect}) max(freqRef{tempSelect})]);
        end
    elseif plotSelect == "Imag"
        for k = 1:nsoc
            ax = axZ(k);
            cla(ax);
            ax.XScale = 'log';
            ax.YScale = 'linear';
            loglog(ax,freqRef{tempSelect},imag(Zref{tempSelect}(:,k)),'b.');
            hold(ax,'on');
            linesZ(k) = semilogx(ax,NaN,NaN,'r-');
            title(ax,sprintf('%.0f%% SOC',socPctRef{tempSelect}(k)));
            ylim(ax,[min(imag(Zref{tempSelect}(:,k))) max(imag(Zref{tempSelect}(:,k)))]);
            xlim(ax,[min(freqRef{tempSelect}) max(freqRef{tempSelect})]);
        end 
    end
    updateUIFig(lastParamValues);
end
end

function ax = formatAxes(ax)
    ax.PlotBoxAspectRatioMode = 'manual';
    ax.PlotBoxAspectRatio = [2/(sqrt(5)-1) 1 1];
    ax.LineWidth = 1;
    ax.FontName = 'Times';
    ax.FontSize = 11;
    ax.TitleFontWeight = 'normal';
    ax.TitleFontSizeMultiplier = 1.2;
    ax.LabelFontSizeMultiplier = 1.1;
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
end