function regressionData = fitLinearEIS(labSpectra,labOCPFit,initialModel)
%FITLINEAREIS Regress linear EIS model to laboratory impedance spectra.
%
% -- Usage --
% regressionData = fitLinearEIS(labSpectra,labOCPFit,initialModel) regresses
%   the reduced-layer Warburg-resistance transfer-function model for LMB 
%   to the linear impedanxe spectra in the structure labSpectra created 
%   with loadLabNLEIS. labOCPFit is the OCP model to use in performing the 
%   EIS regression. initialModel is a model containing initial estimates of 
%   parameter values for the cell.
%
% -- Changelog --
% 2023.06.30 | Created | Wesley Hileman <whileman@uccs.edu>

% Constants.
tplus0 = 0.4;    % Li+ transference number, guess for estimating psi
R = TB.const.R;  % molar gas constant [J/mol/K]
F = TB.const.F;  % Faraday's constant [C/mol]

ocpmodel = MSMR(labOCPFit.MSMR);
ocptest = labOCPFit.ocptest;

% Compute true SOC and lithiation for each SOC setpoint.
zmin = ocpmodel.zmin;
zmax = ocpmodel.zmax;
QdisAhCum = cumsum(labSpectra.QdisAh);
socTrue = 100*(1-QdisAhCum/ocptest.QAh);
theta = zmax + (socTrue/100)*(zmin-zmax);
nsoc = length(socTrue);

% Precompute Uocp, d(Uocp)/d(theta), xj at each SOC setpoint.
% Later used to compute Rct(p) inside optimization loop.
[~,~,~,ocpData] = ocpmodel.ocp('theta',theta,'TdegC',labSpectra.TdegC);

% Convert initial model to reduced-layer Warburg-resistance model.
% Fetch initial values of model parameters.
initialModel = convertCellModel(initialModel,'RLWORM');
initial = getCellParams(initialModel,'TdegC',labSpectra.TdegC);

% Fetch linear impedance measured in the laboratory.
freqLab = labSpectra.lin.freq;
Zlab = labSpectra.lin.Z;
TdegC = labSpectra.TdegC;

% Build model -------------------------------------------------------------
% Known parameters.
params.const.Q = fastopt.param('fix',ocptest.QAh);
params.pos.theta0 = fastopt.param('fix',ocpmodel.zmax);
params.pos.theta100 = fastopt.param('fix',ocpmodel.zmin);
params.pos.X = fastopt.param('fix',ocpmodel.Xj);
params.pos.U0 = fastopt.param('fix',ocpmodel.Uj0);
params.pos.omega = fastopt.param('fix',ocpmodel.Wj);

% Fixed parameters.
% Rdl(p|n) are low for the Sion cells, so fix Rdl(p|n)=0.
params.pos.Rdl = fastopt.param('fix',0);
params.neg.Rdl = fastopt.param('fix',0);
params.neg.Rf = fastopt.param('fix',0);  % Lump into tab resistance.
% Setting wDL=0 results in ideal CPEs for the double-layer capacitors;
% okay since we will not be converting to the time domain inside of the
% optimization routine.
params.pos.wDL = fastopt.param('fix',0);
params.neg.wDL = fastopt.param('fix',0);
% Symmetry factors for the positive electrode are not very identifiable,
% we'll asume alpha=0.5 for all galleries.
params.pos.alpha = fastopt.param('fix',0.5*ones(ocpmodel.J,1));
% Symmetry factor for the negative electrode is not identifiable due to
% linear nature of the small-signal impedance model.
params.neg.alpha = fastopt.param('fix',0.5);

% Electrolyte parameters.
% psi,W,tauW are not separately identifiable, so we fix psi=R/(F*(1-t+0)), 
% the value predicted by the Einstien relationship.
params.const.psi = fastopt.param('fix',R/F/(1-tplus0));
params.const.W = fastopt.param;
params.pos.tauW = fastopt.param('logscale',true);
params.pos.kappa = fastopt.param('logscale',true);
params.eff.tauW = fastopt.param('logscale',true);
params.eff.kappa = fastopt.param('logscale',true);

% Porous electrode parameters.
params.pos.Dsref = fastopt.param('logscale',true);
params.pos.nF = fastopt.param;
params.pos.k0 = fastopt.param('len',ocpmodel.J,'logscale',true);
params.pos.sigma = fastopt.param;
params.pos.Cdl = fastopt.param;
params.pos.nDL = fastopt.param;
params.pos.Rf = fastopt.param;

% Lithium-metal electrode parameters.
params.neg.k0 = fastopt.param('logscale',true);
params.neg.Cdl = fastopt.param;
params.neg.nDL = fastopt.param;

% Cell package parameters.
params.pkg.R0 = fastopt.param;  % Tab resistance.
params.pkg.L0 = fastopt.param;  % Package/cable inductace.

modelspec = fastopt.modelspec(params);


% Define optimization bounds ----------------------------------------------
init = initial;

% Electrolyte parameters.
lb.const.W = 0.1;                   ub.const.W = 10;
lb.eff.tauW = init.eff.tauW/100;    ub.eff.tauW = init.eff.tauW*100;
lb.pos.tauW = init.pos.tauW/100;    ub.pos.tauW = init.pos.tauW*100;
lb.eff.kappa = init.eff.kappa/100;  ub.eff.kappa = init.eff.kappa*100;
lb.pos.kappa = init.pos.kappa/100;  ub.pos.kappa = init.pos.kappa*100;

% Porous-electrode parameters.
lb.pos.Dsref = init.pos.Dsref/1000; ub.pos.Dsref = init.pos.Dsref*1000;
lb.pos.nF = 0.5;                    ub.pos.nF = 1;
lb.pos.k0 = 1e-6*ones(ocpmodel.J,1);ub.pos.k0 = 1e6*ones(ocpmodel.J,1);
lb.pos.sigma = init.pos.sigma/10;   ub.pos.sigma = init.pos.sigma*10;
lb.pos.Cdl = init.pos.Cdl/10;       ub.pos.Cdl = init.pos.Cdl*10;
lb.pos.nDL = 0.5;                   ub.pos.nDL = 1;
lb.pos.Rf = init.pos.Rf/100;        ub.pos.Rf = init.pos.Rf*100;

% Lithium-metal electrode parameters.
lb.neg.k0 = init.neg.k0/1000;       ub.neg.k0 = init.neg.k0*1000;
lb.neg.Cdl = init.neg.Cdl/10;       ub.neg.Cdl = init.neg.Cdl*10;
lb.neg.nDL = 0.5;                   ub.neg.nDL = 1;

% Cell package parameters.
lb.pkg.R0 = 0;                      ub.pkg.R0 = init.pkg.R0*10;
lb.pkg.L0 = 0;                      ub.pkg.L0 = init.pkg.L0*10;

% Let the user adjust the initial values and bounds.
linesNyquist = gobjects(nsoc,1);
lineRctp = gobjects(1,1);
linesRctjp = gobjects(ocpmodel.J,1);
[init,lb,ub] = fastopt.tweekgui( ...
    modelspec,init,lb,ub,@initializeUIPlot,@updateUIPlot);

% Perform regression ------------------------------------------------------

estimate = fastopt.particleswarm( ...
    @cost,modelspec, ...
    fastopt.pack(lb,modelspec,'coerce',true), ...
    fastopt.pack(ub,modelspec,'coerce',true), ...
    'initial',fastopt.pack(init,modelspec,'coerce',true), ...
    'particleCount',5000,'swarmIterations',10, ...
    'fminconIterations',1000,'trackTrajectory',false, ...
    'hybrid','none');
regressionData.estimate = estimate;

[est2,lb2,ub2] = fastopt.tweekgui( ...
    modelspec,estimate,lb,ub,@initializeUIPlot,@updateUIPlot);

function J = cost(model)
    % Calculate impedance predicted by the linear EIS model.
    Zmodel = getModelImpedance(model,freqLab,socTrue,TdegC,ocpData);

    % Compute total residual between model impedance and measured
    % impedance across all spectra.
    J = sum((abs(Zmodel-Zlab)./abs(Zlab)).^2,'all');
end % cost()

function initializeUIPlot(parent,~,~,~)
    gridtop = uigridlayout(parent,[1 2]);
    gridtop.ColumnWidth = {'3x','1x'};
    gridplots = uigridlayout(gridtop,[5 ceil(nsoc/5)]);
    for k = 1:nsoc
        ax = uiaxes(gridplots);
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
        plot(ax,real(Zlab(:,k)),-imag(Zlab(:,k)),'b.');
        hold(ax,'on');
        linesNyquist(k) = plot(ax,NaN,NaN,'r-');
        title(ax,sprintf('%.0f%% SOC',socTrue(k)));
        setAxesNyquist('axes',ax);
    end
    gridsecplots = uigridlayout(gridtop,[2 1]);
    ax = uiaxes(gridsecplots);
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
    lineRctp(1) = semilogy(ax,NaN,NaN,'k-');
    title(ax,'R_{ct}^p vs SOC');
    ax = uiaxes(gridsecplots);
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
    for k = 1:length(linesRctjp)
        linesRctjp(k) = semilogy(ax,NaN,NaN,'-');
        hold(ax,'on');
    end
    title(ax,'R_{ct,j}^p vs SOC');
end % initializeUIPlot()

function updateUIPlot(~,model,~,~)
    % Calculate impedance predicted by the linear EIS model.
    Zmodel = getModelImpedance(model,freqLab,socTrue,TdegC,ocpData);

    % Calculate Rctp predicted by the perturbation model.
    socPct = linspace(0,100,100);
    t = zmax + (socPct/100)*(zmin-zmax);
    pData = getPerturbationResistance( ...
        model,t,'TdegC',TdegC,'ComputeRctj',true);
    Rctp = pData.parts.Rct_p;
    Rctpj = pData.parts.Rctj_p;

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

end % fitLinearEIS()

function Zmodel = getModelImpedance(model,freq,socPct,TdegC,ocpData)
    %GETMODELIMPEDANCE Calculate impedance from TFs at set of frequency and
    %  SOC points.
    nfreq = length(freq);
    nsoc = length(socPct);
    s = 1j*2*pi*freq;  % Laplace variable
    Zmodel = zeros(nfreq,nsoc);

    % Evalulate values of SOC-invariant parameters.
    % Hack to prevent tfLMB from performing expensive calculation
    % of SOC-variant parameters which are already pre-calculated in the
    % ocpData structure.
    tmp = model;
    tmp.pos.Uocp = NaN;
    tmp.pos.dUocp = NaN;
    tfData = tfLMB(s,model,'TdegC',TdegC,'Calc11',false,'Calc22',false);
    baseParam = tfData.param;

    for k = 1:nsoc
        param = baseParam;
        ocpdata = ocpData;
        mod = model;
        % Update values of SOC-variant parameters.
        % - SOC and lithiation: socp, thetap.
        param.socp = socPct(k)/100;
        param.thetap = ocpdata.theta(k);
        % - OCP function: Uocp, dUocp.
        param.Uocpp = ocpdata.Uocp(k);
        param.dUocpp = ocpdata.dUocp(k);
        % - Charge-transfer resistance: Rctp.
        xj = ocpData.xj;
        k0 = mod.pos.k0;
        omega = mod.pos.omega;
        alpha = mod.pos.alpha;
        X = mod.pos.X;
        i0j = k0.*xj(:,k).^(omega.*alpha).*(X-xj(:,k)).^(omega.*(1-alpha));
        i0 = sum(i0j);
        param.Rctp = 1./i0/ocpData.f;
        % - Solid diffusivity: Dsp.
        if isfield(mod.pos,'Dsref')
            param.Dsp = -mod.pos.Dsref*ocpData.f*...
                        param.thetap*(1-param.thetap)*param.dUocpp;
        else
            % Model does not employ SOC-variant solid diffusivity;
            % Ds remains the same at all SOC setpoints.
        end

        % Calculate impedance at this setpoint.
        tfData = tfLMB(s,mod,'ParameterValues',param);
        Zmodel(:,k) = tfData.h11.tfVcell();
    end % for

    % Add impedance contributed by package.
    Zmodel = Zmodel + model.pkg.R0 + model.pkg.L0*s;
end % getModelImpedance()