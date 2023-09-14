function [initializeFcn, updateFcn] = uiRPTCallbacks(modelspec,hcyc,eis)
%UIRPTCALLBACKS

% Constants.
nsoc = length(eis.socPct);
J = modelspec.params.pos__U0.len;

% Globals.
axZ = gobjects(nsoc,1);
gridz = gobjects(1,1);
linesZ = gobjects(nsoc,1);
lastParamValues = [];
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

    % Construct controls.
    gridcontrols = uigridlayout(panelcontrols,[3 1]);
    gridcontrols.RowHeight = {30, 30, '1x'};
    zplotselect = uidropdown(gridcontrols, ...
        "Items",["Nyquist","Bode Magnitude","Bode Phase","Real Part","Imaginary Part"], ...
        "ItemsData",["Nyq","BodeMag","BodePhase","Real","Imag"],...
        "ValueChangedFcn",@(src,event)updateZPlotType(event.Value));

    % Construct impedance plots.
    gridz = uigridlayout(panelz,[ceil(nsoc/5) 5]);
    for k = 1:nsoc
        ax = uiaxes(gridz);
        axZ(k) = ax;
        formatAxes(ax);
        plot(ax,real(eis.lin.Z(:,k)),-imag(eis.lin.Z(:,k)),'b.');
        hold(ax,'on');
        linesZ(k) = plot(ax,NaN,NaN,'r-');
        title(ax,sprintf('%.0f%% SOC',eis.socPct(k)));
        setAxesNyquist('axes',ax);
    end
end % initializeUIPlot()

function updateUIFig(paramValues)
    % Save parameter values so we can update the plots when plot type
    % changes.
    lastParamValues = paramValues;

    % Calculate impedance predicted by the linear EIS model.
    Zmodel = getLinearImpedance( ...
        paramValues,eis.lin.freq,eis.socPct,eis.TdegC);

    % Update model predictions on plots.
    if plotSelect == "Nyq"
        for idxSOC = 1:nsoc
            linesZ(idxSOC).XData = real(Zmodel(:,idxSOC));
            linesZ(idxSOC).YData = -imag(Zmodel(:,idxSOC));
        end
    elseif plotSelect == "BodeMag"
        for idxSOC = 1:nsoc
            linesZ(idxSOC).XData = eis.lin.freq{tempSelect};
            linesZ(idxSOC).YData = abs(Zmodel(:,idxSOC));
        end
    elseif plotSelect == "BodePhase"
        for idxSOC = 1:nsoc
            linesZ(idxSOC).XData = eis.lin.freq{tempSelect};
            linesZ(idxSOC).YData = angle(Zmodel(:,idxSOC))*180/pi;
        end
    elseif plotSelect == "Real"
        for idxSOC = 1:nsoc
            linesZ(idxSOC).XData = eis.lin.freq{tempSelect};
            linesZ(idxSOC).YData = real(Zmodel(:,idxSOC));
        end
    elseif plotSelect == "Imag"
        for idxSOC = 1:nsoc
            linesZ(idxSOC).XData = eis.lin.freq{tempSelect};
            linesZ(idxSOC).YData = imag(Zmodel(:,idxSOC));
        end
    end
end % updateUIPlot()

function updateZPlotType(plottype)
    plotSelect = plottype;
    redrawPlots();
end % updateZPlotType()

function redrawPlots()
    if plotSelect == "Nyq"    
        for k = 1:nsoc
            ax = axZ(k);
            cla(ax);
            ax.XScale = 'linear';
            ax.YScale = 'linear';
            plot(ax,real(eis.lin.Z(:,k)),-imag(eis.lin.Z(:,k)),'b.');
            hold(ax,'on');
            linesZ(k) = plot(ax,NaN,NaN,'r-');
            title(ax,sprintf('%.0f%% SOC',eis.socPct(k)));
            setAxesNyquist('axes',ax);
        end 
    elseif plotSelect == "BodeMag"
        for k = 1:nsoc
            ax = axZ(k);
            cla(ax);
            ax.XScale = 'log';
            ax.YScale = 'log';
            loglog(ax,eis.lin.freq,abs(eis.lin.Z(:,k)),'b.');
            hold(ax,'on');
            linesZ(k) = loglog(ax,NaN,NaN,'r-');
            title(ax,sprintf('%.0f%% SOC',eis.socPct(k)));
            ylim(ax,[min(abs(eis.lin.Z(:,k))) max(abs(eis.lin.Z(:,k)))]);
            xlim(ax,[min(eis.lin.freq) max(eis.lin.freq)]);
        end 
    elseif plotSelect == "BodePhase"
        for k = 1:nsoc
            ax = axZ(k);
            cla(ax);
            ax.XScale = 'log';
            ax.YScale = 'linear';
            semilogx(ax,eis.lin.freq,angle(eis.lin.Z(:,k))*180/pi,'b.');
            hold(ax,'on');
            linesZ(k) = semilogx(ax,NaN,NaN,'r-');
            title(ax,sprintf('%.0f%% SOC',eis.socPct(k)));
            ylim(ax,[min(angle(eis.lin.Z(:,k))*180/pi) max(angle(eis.lin.Z(:,k))*180/pi)]);
            xlim(ax,[min(eis.lin.freq) max(eis.lin.freq)]);
        end
    elseif plotSelect == "Real"
        for k = 1:nsoc
            ax = axZ(k);
            cla(ax);
            ax.XScale = 'log';
            ax.YScale = 'linear';
            loglog(ax,eis.lin.freq,real(eis.lin.Z(:,k)),'b.');
            hold(ax,'on');
            linesZ(k) = semilogx(ax,NaN,NaN,'r-');
            title(ax,sprintf('%.0f%% SOC',eis.socPct(k)));
            ylim(ax,[min(real(eis.lin.Z(:,k))) max(real(eis.lin.Z(:,k)))]);
            xlim(ax,[min(eis.lin.freq) max(eis.lin.freq)]);
        end
    elseif plotSelect == "Imag"
        for k = 1:nsoc
            ax = axZ(k);
            cla(ax);
            ax.XScale = 'log';
            ax.YScale = 'linear';
            loglog(ax,eis.lin.freq,imag(eis.lin.Z(:,k)),'b.');
            hold(ax,'on');
            linesZ(k) = semilogx(ax,NaN,NaN,'r-');
            title(ax,sprintf('%.0f%% SOC',eis.socPct(k)));
            ylim(ax,[min(imag(eis.lin.Z(:,k))) max(imag(eis.lin.Z(:,k)))]);
            xlim(ax,[min(eis.lin.freq) max(eis.lin.freq)]);
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
