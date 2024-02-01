function data = fitTF(regspec,varargin)
%FITTF Regress transfer-function model to electrochemical impedace spectra.
%
% -- Usage --
% data = FITTF(regspec)
%
% -- Changelog --
% 2024.01.29 | Cleanup, enable use of synthetic data | WH
% 2023.09.19 | Fix error computing theta, wrong theta min/max | Wes H.
% 2023.08.27 | Store separate IS for Ds/i0 based on distance | Wes H.
% 2023.08.10 | Add option to save intermediate solutions (IS) | Wes H.
% 2023.07.26 | Update for multiple temperatures | Wes H.
% 2023.07.24 | Update for multiple diffusion and kinetics models | Wes H.
% 2023.06.30 | Created | Wesley Hileman <whileman@uccs.edu>

%FITLINEAREIS Regress linear EIS model to laboratory impedance spectra.
%
% -- Usage --
% regressionData = fitLinearEIS(labSpectra,labOCPFit,initialModel) regresses
%   the reduced-layer Warburg-resistance transfer-function model for LMB 
%   to the linear impedance spectra in the structure labSpectra created 
%   with loadLabNLEIS. labOCPFit is the OCP model to use in performing the 
%   EIS regression. initialModel is a model containing initial estimates of 
%   parameter values for the cell.
%
% Supply datasets at multiple temperatures by specifying LABSPECTRA as a 
% structure array, containing a scalar structure for each temperature 
% dataset.
%
% Output structure:
%   .estimate    Cell model containing estimated parameter values
%   .modelspec   fastopt specification for EIS model
%   .values      Structure of estimated parameter values
%   .initial     Structure of initial parameter values
%   .lb          Structure of lower bounds for parameters
%   .ub          Structure of upper bounds for parameters
%   .Zmodel      Matrix of model-predicted impedance: d1=frequency d2=soc
%   .Zlab        Matrix of lab impedance: d1=frequency d2=soc
%   .freq        Cyclic frequency vector
%   .socPctTrue  SOC vector, calculated from QdisAh
%   .TdegC       Temperature
%   .arg         Structure of arguments supplied to the function.
%

parser = inputParser;
parser.addRequired('regspec',@(x)isstruct(x)&&isscalar(x));
parser.addParameter('WeightFcn',[],@(x)isa(x,'function_handle'));
parser.addParameter('UseParallel',true,@(x)islogical(x)&&isscalar(x));
% --
% NOTE: IntermediateSoln stuff is disabled when UseParallel option is true;
% MATLAB can't easily share data-structures between threads.
parser.addParameter('IntermediateSolns',100,@(x)isnumeric(x)&&isscalar(x)); % num previous minimizers to keep
parser.addParameter('IntermediateSolnEpsilon',struct('Ds',0.03,'k0',0.03),@(x)isscalar(x)&&istruct(x));
% --
parser.parse(regspec,varargin{:});
arg = parser.Results;  % structure of validated arguments


% Collect arguments -------------------------------------------------------

spec = arg.regspec;
TdegC = spec.TdegC;  % temperature vector [degC].
ntemp = spec.ntemp;
eis = spec.eis;
ocp = spec.ocp;
init = spec.init;
lb = spec.lb;
ub = spec.ub;
modelspec = spec.modelspec;

% Collect linear impedance measured in the laboratory.
clear lab;
for m = ntemp:-1:1
    tmp = eis(m);
    tmp.theta = ocp.theta0 + (eis(m).socPct/100)*(ocp.theta100 - ocp.theta0);
    tmp.ocpData = MSMR(ocp).ocp('theta',tmp.theta,'TdegC',TdegC(m));
    lab(m) = tmp;
end 

% Collect weight matrix.
weights = cell(ntemp,1);
for m = 1:ntemp
    weights{m} = ones(length(lab(m).freq),length(lab(m).socPct));
    if ~isempty(arg.WeightFcn)
        for idxSOC = 1:length(lab(m).socPct)
            for idxFreq = 1:length(lab(m).freq)
                weights{m}(idxFreq,idxSOC) = arg.WeightFcn( ...
                    lab(m).freq(idxFreq),lab(m).socPct(idxSOC),TdegC(m));
            end % for
        end % for
    end % if
end % for


% Perform regression ------------------------------------------------------

% Allocate storage for intermediate solutions.
topJ.Ds = inf(1,arg.IntermediateSolns);
topJ.k0 = inf(1,arg.IntermediateSolns);
topSoln.Ds(arg.IntermediateSolns) = init;
topSoln.k0(arg.IntermediateSolns) = init;

% Perform regression.
[plotInit, plotUpdate] = uiImpedanceCallbacks(modelspec,lab);
data = fastopt.uiparticleswarm(@cost,modelspec,init,lb,ub, ...
    'PlotInitializeFcn',plotInit,'PlotUpdateFcn',plotUpdate, ...
    'UseParallel',arg.UseParallel);

% Collect output data.
data.topJ = topJ;
data.topSoln = topSoln;
data.modelspec = modelspec;
models = fastopt.splittemps(data.values,modelspec);
clear fit;
for m = ntemp:-1:1
    tmp.freq = lab(m).freq;
    tmp.socPct = lab(m).socPct;
    tmp.TdegC = lab(m).TdegC;
    tmp.Z = calcZ(models(m),lab(m));
    fit(m) = tmp;
end
data.fit = fit;
data.lab = lab;
data.TdegC = TdegC;
data.argpso__ = data.arg;
data = rmfield(data,{'origin__','arg'});
data.arg = arg;
data.origin__ = 'fitTF';

function J = cost(model)
    modelVect = fastopt.splittemps(model,modelspec);

    % Ignore singular matrix warnings in solving TF model (usu. occurs at high
    % frequencies). TODO: find high-frequency TF solution.
    warning('off','MATLAB:nearlySingularMatrix');

    % Accumulate cost at each temperature.
    J = 0;
    for km = 1:ntemp
        % Calculate impedance predicted by the linear EIS model.
        Zmodel = calcZ(modelVect(km),lab(km));
        % Compute total residual between model impedance and measured
        % impedance across all spectra.
        J = J + sum(weights{km}.*(abs(Zmodel-lab(km).Z)./abs(lab(km).Z)).^2,'all');
    end

    % Re-enable singular matrix warning.
    warning('on','MATLAB:nearlySingularMatrix');

    % Store intermediate solutions.
    if ~arg.UseParallel
        Jlt = J<topJ.Ds; 
        if any(Jlt) && strcmpi(spec.type,'lut')
            if isinf(topJ.Ds(1))
                store = true;
            else
                pos = [topSoln.Ds(~isinf(topJ.Ds)).pos];
                n0Ds = lognormalize( ...
                    [pos.DsLinear],lb.pos.DsLinear,ub.pos.DsLinear);
                nDs = lognormalize( ...
                    model.pos.DsLinear,lb.pos.DsLinear,ub.pos.DsLinear);
                store = all(vecnorm(nDs-n0Ds)>arg.IntermediateSolnEpsilon.Ds);
            end % else
            if store
                ind = find(Jlt,1,'first');
                topJ.Ds = [topJ.Ds(1:ind-1) J topJ.Ds(ind:end-1)];
                topSoln.Ds = [topSoln.Ds(1:ind-1) model topSoln.Ds(ind:end-1)];
            end % if
        end % if
        Jlt = J<topJ.k0; 
        if any(Jlt) && strcmpi(spec.type,'lut')
            if isinf(topJ.k0(1))
                store = true;
            else
                pos = [topSoln.k0(~isinf(topJ.k0)).pos];
                n0k0 = lognormalize( ...
                    [pos.k0Linear],lb.pos.k0Linear,ub.pos.k0Linear);
                nk0 = lognormalize( ...
                    model.pos.k0Linear,lb.pos.k0Linear,ub.pos.k0Linear);
                store = all(vecnorm(nk0-n0k0)>arg.IntermediateSolnEpsilon.k0);
            end % else
            if store
                ind = find(Jlt,1,'first');
                topJ.k0 = [topJ.k0(1:ind-1) J topJ.k0(ind:end-1)];
                topSoln.k0 = [topSoln.k0(1:ind-1) model topSoln.k0(ind:end-1)];
            end % if
        end % if
    end % if
end % cost()

function n = lognormalize(x,lb,ub)
    n = (log10(x)-log10(lb))./(log10(ub)-log10(lb));
end

end % fitLinearEIS()


% Utility Functions -------------------------------------------------------

function Z = calcZ(model,data)
%CALCZ Calculate TF impedance at given frequency / SOC points.
%
% -- Usage --
% Z = getLinearImpedance(cellParams,data)
%
% -- Input --
% cellParams = struct of cell parameter values (fast) OR cell model struct 
%   which will be evalulated at specified temperature (slow)
%
% data = struct with the following fields:
%   .freq    : frequency vector [Hz]
%   .socPct  : SOC vector [%]
%   .TdegC   : temperature [degC]
%   .ocpData : output of MSMR.ocp() @ TdegC and socPct vector 
%              (OPTIONAL, reduces computation time)
%
% -- Output --
% Z = matrix of complex impedance values (dim1<=>freq, dim2<=>socPct) [Ohm]

% freq,socPct,TdegC,ocpData

freq = data.freq;
socPct = data.socPct;
TdegC = data.TdegC;

nfreq = length(freq);
nsoc = length(socPct);
s = 1j*2*pi*freq;  % Laplace variable

if isCellModel(model)
    params = getCellParams(model,'TdegC',TdegC);
else
    params = model;
end

useCachedOCP = false;
if isfield(data,'ocpData') && ~isempty(data.ocpData)
    ocpModel = MSMR(data.ocpData);
    ocpData = data.ocpData;
    ctData = ocpModel.RctCachedOCP(params.pos,data.ocpData);
    dsData = ocpModel.DsCachedOCP(params.pos,data.ocpData);
    useCachedOCP = true;
end

% Calculate impedance at each SOC setpoint.
Z = zeros(nfreq,nsoc);
for k = 1:nsoc
    if useCachedOCP
        % Update OCP, charge-transfer resistance, and spolid diffusivity 
        % at each SOC setpoint  manually for improved performance.
        params.pos.Uocp = ocpData.Uocp(k);
        params.pos.dUocp = ocpData.dUocp(k);
        params.pos.Rct = ctData.Rct(k);
        params.pos.Ds = dsData.Ds(k);
    end
    tfData = tfLMB(s,params,'TdegC',TdegC,'socPct',socPct(k));
    Z(:,k) = tfData.h11.tfZcell();
end % for

% Add impedance contributed by package.
if isfield(params,'pkg')
    Z = Z + params.pkg.R0 + params.pkg.L0*s;
end
end


% UI Utility Functions ----------------------------------------------------

function [initializeFcn, updateFcn] = uiImpedanceCallbacks(modelspec,ref)
%UIIMPEDANCECALLBACKS

%Zref,freqRef,socPctRef,TdegC
TdegC = [ref.TdegC];

nsoc = length(ref(1).socPct);
J = modelspec.params.pos__U0.len;

% Define globals.
axZ = gobjects(nsoc,1);
zplotselect = gobjects(1,1);
tempselect = gobjects(1,1);
gridz = gobjects(1,1);
linesZ = gobjects(nsoc,1);
lineDsp = gobjects(1,1);
lineDspInterp = gobjects(1,1);
linei0p = gobjects(1,1);
linei0pInterp = gobjects(1,1);
lineRctp = gobjects(1,1);
linesRctjp = gobjects(J,1);
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
        plot(ax,real(ref(tempSelect).Z(:,k)),-imag(ref(tempSelect).Z(:,k)),'b.');
        hold(ax,'on');
        linesZ(k) = plot(ax,NaN,NaN,'r-');
        title(ax,sprintf('%.0f%% SOC',ref(tempSelect).socPct(k)));
        setAxesNyquist('axes',ax);
    end

    % Construct charge-transfer resistance plots.
    gridauxplots = uigridlayout(panelct,[4 1]);
    ax = uiaxes(gridauxplots);
    formatAxes(ax);
    lineDsp(1) = semilogy(ax,NaN,NaN,'k-');
    hold(ax,'on');
    lineDspInterp(1) = semilogy(ax,NaN,NaN,'ro');
    title(ax,'D_{s}^p vs SOC')
    ax = uiaxes(gridauxplots);
    formatAxes(ax);
    linei0p(1) = semilogy(ax,NaN,NaN,'k-');
    hold(ax,'on');
    linei0pInterp(1) = semilogy(ax,NaN,NaN,'ro');
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

    % Select appropriate data-set.
    paramValues = fastopt.splittemps(paramValues,modelspec);
    paramValues = paramValues(tempSelect);

    % Calculate impedance predicted by the linear EIS model.
    Zmodel = calcZ(paramValues,ref(tempSelect));

    % Calculate Rctp and Dsp.
    socPct = linspace(0,100,100);
    theta0 = paramValues.pos.theta0;
    theta100 = paramValues.pos.theta100;
    t = theta0 + (socPct/100)*(theta100-theta0);
    if isfield(paramValues.pos,'DsTheta')
        DsTheta = paramValues.pos.DsTheta;
        DsSOCPct = 100*(DsTheta-theta0)/(theta100-theta0);
    end
    if isfield(paramValues.pos,'k0Theta')
        k0Theta = paramValues.pos.k0Theta;
        k0SOCPct = 100*(k0Theta-theta0)/(theta100-theta0);
    end
    ocpmodel = MSMR(paramValues.pos);
    ctData = ocpmodel.Rct(paramValues.pos,'theta',t,'TdegC',TdegC(tempSelect));
    dsData = ocpmodel.Ds(paramValues.pos,'theta',t,'TdegC',TdegC(tempSelect));
    Rctp = ctData.Rct;
    Rctpj = ctData.Rctj;
    i0 = ctData.i0;
    Ds = dsData.Ds;

    % Update model predictions on plots.
    if plotSelect == "Nyq"
        for idxSOC = 1:nsoc
            linesZ(idxSOC).XData = real(Zmodel(:,idxSOC));
            linesZ(idxSOC).YData = -imag(Zmodel(:,idxSOC));
        end
    elseif plotSelect == "BodeMag"
        for idxSOC = 1:nsoc
            linesZ(idxSOC).XData = ref(tempSelect).freq;
            linesZ(idxSOC).YData = abs(Zmodel(:,idxSOC));
        end
    elseif plotSelect == "BodePhase"
        for idxSOC = 1:nsoc
            linesZ(idxSOC).XData = ref(tempSelect).freq;
            linesZ(idxSOC).YData = angle(Zmodel(:,idxSOC))*180/pi;
        end
    elseif plotSelect == "Real"
        for idxSOC = 1:nsoc
            linesZ(idxSOC).XData = ref(tempSelect).freq;
            linesZ(idxSOC).YData = real(Zmodel(:,idxSOC));
        end
    elseif plotSelect == "Imag"
        for idxSOC = 1:nsoc
            linesZ(idxSOC).XData = ref(tempSelect).freq;
            linesZ(idxSOC).YData = imag(Zmodel(:,idxSOC));
        end
    end
    lineDsp.XData = socPct;
    lineDsp.YData = Ds;
    if isfield(paramValues.pos,'DsLinear')
        lineDspInterp.XData = DsSOCPct;
        lineDspInterp.YData = paramValues.pos.DsLinear;
    elseif isfield(paramValues.pos,'DsSpline')
        lineDspInterp.XData = DsSOCPct;
        lineDspInterp.YData = paramValues.pos.DsSpline;
    end
    lineRctp.XData = socPct;
    lineRctp.YData = Rctp;
    linei0p.XData = socPct;
    linei0p.YData = i0;
    if isfield(paramValues.pos,'k0Linear')
        linei0pInterp.XData = k0SOCPct;
        linei0pInterp.YData = paramValues.pos.k0Linear;
    elseif isfield(paramValues.pos,'k0Spline')
        linei0pInterp.XData = k0SOCPct;
        linei0pInterp.YData = paramValues.pos.k0Spline;
    end
    for j = 1:J
        linesRctjp(j).XData = socPct;
        linesRctjp(j).YData = Rctpj(j,:);
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
            plot(ax,real(ref(tempSelect).Z(:,k)),-imag(ref(tempSelect).Z(:,k)),'b.');
            hold(ax,'on');
            linesZ(k) = plot(ax,NaN,NaN,'r-');
            title(ax,sprintf('%.0f%% SOC',ref(tempSelect).socPct(k)));
            setAxesNyquist('axes',ax);
        end 
    elseif plotSelect == "BodeMag"
        for k = 1:nsoc
            ax = axZ(k);
            cla(ax);
            ax.XScale = 'log';
            ax.YScale = 'log';
            loglog(ax,ref(tempSelect).freq,abs(ref(tempSelect).Z(:,k)),'b.');
            hold(ax,'on');
            linesZ(k) = loglog(ax,NaN,NaN,'r-');
            title(ax,sprintf('%.0f%% SOC',ref(tempSelect).socPct(k)));
            ylim(ax,[min(abs(ref(tempSelect).Z(:,k))) max(abs(ref(tempSelect).Z(:,k)))]);
            xlim(ax,[min(ref(tempSelect).freq) max(ref(tempSelect).freq)]);
        end 
    elseif plotSelect == "BodePhase"
        for k = 1:nsoc
            ax = axZ(k);
            cla(ax);
            ax.XScale = 'log';
            ax.YScale = 'linear';
            semilogx(ax,ref(tempSelect).freq,angle(ref(tempSelect).Z(:,k))*180/pi,'b.');
            hold(ax,'on');
            linesZ(k) = semilogx(ax,NaN,NaN,'r-');
            title(ax,sprintf('%.0f%% SOC',ref(tempSelect).socPct(k)));
            ylim(ax,[min(angle(ref(tempSelect).Z(:,k))*180/pi) max(angle(ref(tempSelect).Z(:,k))*180/pi)]);
            xlim(ax,[min(ref(tempSelect).freq) max(ref(tempSelect).freq)]);
        end
    elseif plotSelect == "Real"
        for k = 1:nsoc
            ax = axZ(k);
            cla(ax);
            ax.XScale = 'log';
            ax.YScale = 'linear';
            loglog(ax,ref(tempSelect).freq,real(ref(tempSelect).Z(:,k)),'b.');
            hold(ax,'on');
            linesZ(k) = semilogx(ax,NaN,NaN,'r-');
            title(ax,sprintf('%.0f%% SOC',ref(tempSelect).socPct(k)));
            ylim(ax,[min(real(ref(tempSelect).Z(:,k))) max(real(ref(tempSelect).Z(:,k)))]);
            xlim(ax,[min(ref(tempSelect).freq) max(ref(tempSelect).freq)]);
        end
    elseif plotSelect == "Imag"
        for k = 1:nsoc
            ax = axZ(k);
            cla(ax);
            ax.XScale = 'log';
            ax.YScale = 'linear';
            loglog(ax,ref(tempSelect).freq,imag(ref(tempSelect).Z(:,k)),'b.');
            hold(ax,'on');
            linesZ(k) = semilogx(ax,NaN,NaN,'r-');
            title(ax,sprintf('%.0f%% SOC',ref(tempSelect).socPct(k)));
            ylim(ax,[min(imag(ref(tempSelect).Z(:,k))) max(imag(ref(tempSelect).Z(:,k)))]);
            xlim(ax,[min(ref(tempSelect).freq) max(ref(tempSelect).freq)]);
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