function [values,lb,ub] = uiAdjustLinearImpedance( ...
    modelspec,values,lb,ub,Zref,freqRef,socPctRef,TdegC)
%UIADJUSTLINEARIMPEDANCE Launch GUI for interactively adjusting parameters
%  of linear impedace model.
%
% -- Changelog --
% 2023.07.02 | Created | Wesley Hileman <whileman@uccs.edu>

nsoc = length(socPctRef);
J = length(values.pos.U0);
linesNyquist = gobjects(nsoc,1);
lineRctp = gobjects(1,1);
linesRctjp = gobjects(J,1);

[values,lb,ub] = fastopt.tweekgui( ...
    modelspec,values,lb,ub,@initializeUIPlot,@updateUIPlot);

    function initializeUIPlot(parent,~,~,~)
        gridtop = uigridlayout(parent,[1 2]);
        gridtop.ColumnWidth = {'3x','1x'};
        gridplots = uigridlayout(gridtop,[5 ceil(nsoc/5)]);
        for k = 1:nsoc
            ax = uiaxes(gridplots);
            formatAxes(ax);
            plot(ax,real(Zref(:,k)),-imag(Zref(:,k)),'b.');
            hold(ax,'on');
            linesNyquist(k) = plot(ax,NaN,NaN,'r-');
            title(ax,sprintf('%.0f%% SOC',socPctRef(k)));
            setAxesNyquist('axes',ax);
        end
        gridsecplots = uigridlayout(gridtop,[2 1]);
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
    
    function updateUIPlot(~,paramValues,~,~)
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
        for idxSOC = 1:nsoc
            linesNyquist(idxSOC).XData = real(Zmodel(:,idxSOC));
            linesNyquist(idxSOC).YData = -imag(Zmodel(:,idxSOC));
        end
        lineRctp.XData = socPct;
        lineRctp.YData = Rctp;
        for j = 1:ocpmodel.J
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