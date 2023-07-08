function [initializeFcn, updateFcn] = uiImpedanceCallbacks(modelspec,Zref,freqRef,socPctRef,TdegC)
%UIIMPEDANCECALLBACKS

nsoc = length(socPctRef);
J = modelspec.params.pos__U0.len;

% Define globals.
axZ = gobjects(nsoc,1);
zplotselect = gobjects(1,1);
gridz = gobjects(1,1);
linesZ = gobjects(nsoc,1);
lineDsp = gobjects(1,1);
lineDspSpline = gobjects(1,1);
linei0p = gobjects(1,1);
linei0pSpline = gobjects(1,1);
lineRctp = gobjects(1,1);
linesRctjp = gobjects(J,1);
lastParamValues = [];

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
    labinfo = uilabel(gridcontrols,"Text",sprintf('T=%.0fdegC',TdegC));

    % Construct impedance plots.
    gridz = uigridlayout(panelz,[5 ceil(nsoc/5)]);
    for k = 1:nsoc
        ax = uiaxes(gridz);
        axZ(k) = ax;
        formatAxes(ax);
        plot(ax,real(Zref(:,k)),-imag(Zref(:,k)),'b.');
        hold(ax,'on');
        linesZ(k) = plot(ax,NaN,NaN,'r-');
        title(ax,sprintf('%.0f%% SOC',socPctRef(k)));
        setAxesNyquist('axes',ax);
    end

    % Construct charge-transfer resistance plots.
    gridauxplots = uigridlayout(panelct,[4 1]);
    ax = uiaxes(gridauxplots);
    formatAxes(ax);
    lineDsp(1) = semilogy(ax,NaN,NaN,'k-');
    hold(ax,'on');
    lineDspSpline(1) = semilogy(ax,NaN,NaN,'ro');
    title(ax,'D_{s}^p vs SOC')
    ax = uiaxes(gridauxplots);
    formatAxes(ax);
    linei0p(1) = semilogy(ax,NaN,NaN,'k-');
    hold(ax,'on');
    linei0pSpline(1) = semilogy(ax,NaN,NaN,'ro');
    title(ax,'i_{0}^p vs SOC');
    ax = uiaxes(gridauxplots);
    formatAxes(ax);
    lineRctp(1) = semilogy(ax,NaN,NaN,'k-');
    title(ax,'R_{ct}^p vs SOC');
    ax = uiaxes(gridauxplots);
    formatAxes(ax);
    for k = 1:length(linesRctjp)
        linesRctjp(k) = semilogy(ax,NaN,NaN,'-');
        hold(ax,'on');
    end
    title(ax,'R_{ct,j}^p vs SOC');
end % initializeUIPlot()

function updateUIFig(paramValues)
    % Save parameter values so we can update the plots when plot type
    % changes.
    lastParamValues = paramValues;

    % Calculate impedance predicted by the linear EIS model.
    Zmodel = getLinearImpedance(paramValues,freqRef,socPctRef,TdegC);

    % Calculate Rctp and Dsp.
    DsSplineTheta = paramValues.pos.DsSplineTheta;
    k0SplineTheta = paramValues.pos.k0SplineTheta;
    socPct = linspace(0,100,100);
    theta0 = paramValues.pos.theta0;
    theta100 = paramValues.pos.theta100;
    t = theta0 + (socPct/100)*(theta100-theta0);
    DsSplineSOCPct = 100*(DsSplineTheta-theta0)/(theta100-theta0);
    k0SplineSOCPct = 100*(k0SplineTheta-theta0)/(theta100-theta0);
    ocpmodel = MSMR(paramValues.pos);
    ctData = ocpmodel.Rct(paramValues.pos,'theta',t,'TdegC',TdegC);
    dsData = ocpmodel.Ds(paramValues.pos,'theta',t,'TdegC',TdegC);
    Rctp = ctData.Rct;
    Rctpj = ctData.Rctj;
    i0 = ctData.i0;
    Ds = dsData.Ds;

    % Update model predictions on plots.
    plottype = zplotselect.Value;
    if plottype == "Nyq"
        for idxSOC = 1:nsoc
            linesZ(idxSOC).XData = real(Zmodel(:,idxSOC));
            linesZ(idxSOC).YData = -imag(Zmodel(:,idxSOC));
        end
    elseif plottype == "BodeMag"
        for idxSOC = 1:nsoc
            linesZ(idxSOC).XData = freqRef;
            linesZ(idxSOC).YData = abs(Zmodel(:,idxSOC));
        end
    elseif plottype == "BodePhase"
        for idxSOC = 1:nsoc
            linesZ(idxSOC).XData = freqRef;
            linesZ(idxSOC).YData = angle(Zmodel(:,idxSOC))*180/pi;
        end
    elseif plottype == "Real"
        for idxSOC = 1:nsoc
            linesZ(idxSOC).XData = freqRef;
            linesZ(idxSOC).YData = real(Zmodel(:,idxSOC));
        end
    elseif plottype == "Imag"
        for idxSOC = 1:nsoc
            linesZ(idxSOC).XData = freqRef;
            linesZ(idxSOC).YData = imag(Zmodel(:,idxSOC));
        end
    end
    lineDsp.XData = socPct;
    lineDsp.YData = Ds;
    lineDspSpline.XData = DsSplineSOCPct;
    lineDspSpline.YData = paramValues.pos.DsSpline;
    lineRctp.XData = socPct;
    lineRctp.YData = Rctp;
    linei0p.XData = socPct;
    linei0p.YData = i0;
    linei0pSpline.XData = k0SplineSOCPct;
    linei0pSpline.YData = paramValues.pos.k0Spline;
    for j = 1:J
        linesRctjp(j).XData = socPct;
        linesRctjp(j).YData = Rctpj(j,:);
    end
end % updateUIPlot()

function updateZPlotType(plottype)
    if plottype == "Nyq"    
        for k = 1:nsoc
            ax = axZ(k);
            cla(ax);
            ax.XScale = 'linear';
            ax.YScale = 'linear';
            plot(ax,real(Zref(:,k)),-imag(Zref(:,k)),'b.');
            hold(ax,'on');
            linesZ(k) = plot(ax,NaN,NaN,'r-');
            title(ax,sprintf('%.0f%% SOC',socPctRef(k)));
            setAxesNyquist('axes',ax);
        end 
    elseif plottype == "BodeMag"
        for k = 1:nsoc
            ax = axZ(k);
            cla(ax);
            ax.XScale = 'log';
            ax.YScale = 'log';
            loglog(ax,freqRef,abs(Zref(:,k)),'b.');
            hold(ax,'on');
            linesZ(k) = loglog(ax,NaN,NaN,'r-');
            title(ax,sprintf('%.0f%% SOC',socPctRef(k)));
            ylim(ax,[min(abs(Zref(:,k))) max(abs(Zref(:,k)))]);
            xlim(ax,[min(freqRef) max(freqRef)]);
        end 
    elseif plottype == "BodePhase"
        for k = 1:nsoc
            ax = axZ(k);
            cla(ax);
            ax.XScale = 'log';
            ax.YScale = 'linear';
            semilogx(ax,freqRef,angle(Zref(:,k))*180/pi,'b.');
            hold(ax,'on');
            linesZ(k) = semilogx(ax,NaN,NaN,'r-');
            title(ax,sprintf('%.0f%% SOC',socPctRef(k)));
            ylim(ax,[min(angle(Zref(:,k))*180/pi) max(angle(Zref(:,k))*180/pi)]);
            xlim(ax,[min(freqRef) max(freqRef)]);
        end
    elseif plottype == "Real"
        for k = 1:nsoc
            ax = axZ(k);
            cla(ax);
            ax.XScale = 'log';
            ax.YScale = 'linear';
            loglog(ax,freqRef,real(Zref(:,k)),'b.');
            hold(ax,'on');
            linesZ(k) = semilogx(ax,NaN,NaN,'r-');
            title(ax,sprintf('%.0f%% SOC',socPctRef(k)));
            ylim(ax,[min(real(Zref(:,k))) max(real(Zref(:,k)))]);
            xlim(ax,[min(freqRef) max(freqRef)]);
        end
    elseif plottype == "Imag"
        for k = 1:nsoc
            ax = axZ(k);
            cla(ax);
            ax.XScale = 'log';
            ax.YScale = 'linear';
            loglog(ax,freqRef,imag(Zref(:,k)),'b.');
            hold(ax,'on');
            linesZ(k) = semilogx(ax,NaN,NaN,'r-');
            title(ax,sprintf('%.0f%% SOC',socPctRef(k)));
            ylim(ax,[min(imag(Zref(:,k))) max(imag(Zref(:,k)))]);
            xlim(ax,[min(freqRef) max(freqRef)]);
        end 
    end
    updateUIFig(lastParamValues);
end % updateZPlotType()
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