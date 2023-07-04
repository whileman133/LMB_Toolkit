function [initializeFcn, updateFcn] = uiImpedanceCallbacks(modelspec,Zref,freqRef,socPctRef,TdegC)
%UIIMPEDANCECALLBACKS

% Define globals.
nsoc = length(socPctRef);
J = modelspec.params.pos__U0.len;
axZ = gobjects(nsoc,1);
zplotselect = gobjects(1,1);
gridz = gobjects(1,1);
linesZ = gobjects(nsoc,1);
lineRctp = gobjects(1,1);
linesRctjp = gobjects(J,1);
lastParamValues = [];

initializeFcn = @initializeUIFig;
updateFcn = @updateUIFig;

function initializeUIFig(parent)
    gridtop = uigridlayout(parent,[1 2]);
    gridtop.ColumnWidth = {'3x','1x'};
    gridz = uigridlayout(gridtop,[5 ceil(nsoc/5)]);
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
    gridsec = uigridlayout(gridtop,[2 1]);
    gridsec.RowHeight = {50, '1x'};
    zplotselect = uidropdown(gridsec, ...
        "Items",["Nyquist","Bode Magnitude","Bode Phase","Real Part","Imaginary Part"], ...
        "ItemsData",["Nyq","BodeMag","BodePhase","Real","Imag"],...
        "ValueChangedFcn",@(src,event)updateZPlotType(event.Value));
    gridsecplots = uigridlayout(gridsec,[2 1]);
    ax = uiaxes(gridsecplots);
    formatAxes(ax);
    lineRctp(1) = semilogy(ax,NaN,NaN,'k-');
    title(ax,'R_{ct}^p vs SOC');
    ax = uiaxes(gridsecplots);
    formatAxes(ax);
    for k = 1:length(linesRctjp)
        linesRctjp(k) = semilogy(ax,NaN,NaN,'-');
        hold(ax,'on');
    end
    title(ax,'R_{ct,j}^p vs SOC');
end % initializeUIPlot()

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
end

function updateUIFig(paramValues)
    % Save parameter values so we can update the plots when plot type
    % changes.
    lastParamValues = paramValues;

    % Calculate impedance predicted by the linear EIS model.
    Zmodel = getLinearImpedance(paramValues,freqRef,socPctRef,TdegC);

    % Calculate Rctp.
    socPct = linspace(0,100,100);
    theta0 = paramValues.pos.theta0;
    theta100 = paramValues.pos.theta100;
    t = theta0 + (socPct/100)*(theta100-theta0);
    ocpmodel = MSMR(paramValues.pos);
    ctData = ocpmodel.Rct(paramValues.pos,'theta',t,'TdegC',TdegC);
    Rctp = ctData.Rct;
    Rctpj = ctData.Rctj;

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
    lineRctp.XData = socPct;
    lineRctp.YData = Rctp;
    for j = 1:J
        linesRctjp(j).XData = socPct;
        linesRctjp(j).YData = Rctpj(j,:);
    end
end % updateUIPlot()
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