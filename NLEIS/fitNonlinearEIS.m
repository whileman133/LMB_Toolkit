function data = fitNonlinearEIS(labSpectra,labLinearFit,varargin)
%FITNONLINEAREIS Regress nonlinear EIS model to lab impedance spectra.
%
% -- Usage --
% regressionData = fitNonlinearEIS(labSpectra,labLinearFit) regresses
%   the reduced-layer Warburg-resistance second-harmonic impedance model 
%   for LMB to the second-harmonoc impedance spectra in the structure
%   labSpectra created with loadLabNLEIS. labLinearFit is the data
%   structure from the linear EIS regression.
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
% -- Changelog --
% 2023.07.17 | Created | Wesley Hileman <whileman@uccs.edu>

parser = inputParser;
parser.addRequired('labSpectra',@isstruct);
parser.addRequired('labLinearFit',@isstruct);
parser.addParameter('WeightFcn',[],@(x)isa(x,'function_handle'));
parser.addParameter('UseParallel',true);
parser.parse(labSpectra,labLinearFit,varargin{:});
arg = parser.Results;  % structure of validated arguments

% Determine the multiplicity of the dataset, i.e. the number of
% temperatures included.
multiplicity = length(arg.labSpectra);

% Collect NL impedance measured in the laboratory.
freqLab = cell(multiplicity,1);
Zlab = cell(multiplicity,1);
for m = 1:multiplicity
    freqLab{m} = arg.labSpectra(m).h2.freq;
    Zlab{m} = arg.labSpectra(m).h2.Z;
end
TdegC = [arg.labSpectra.TdegC].';  % temperature vector [degC].

% Fetch OCP parameters fit to laboratory data.
ocpmodel = MSMR(labLinearFit.values.pos);
zmin = ocpmodel.zmin;
zmax = ocpmodel.zmax;
socPctTrue = labLinearFit.socPctTrue;
thetaTrue = socPctTrue;
for k = 1:length(socPctTrue)
    thetaTrue{k} = zmax + (socPctTrue{k}/100)*(zmin-zmax);
end
ocpData = cell(multiplicity,1);
for m = 1:multiplicity
    ocpData{m} = ocpmodel.ocp('theta',thetaTrue{m},'TdegC',TdegC(m));
end

% Build model -------------------------------------------------------------
v = labLinearFit.values;

% Known parameters.
% OCP.
params.const.Q = fastopt.param('fix',v.const.Q);
params.pos.theta0 = fastopt.param('fix',v.pos.theta0);
params.pos.theta100 = fastopt.param('fix',v.pos.theta100);
params.pos.X = fastopt.param('fix',v.pos.X);
params.pos.U0 = fastopt.param('fix',v.pos.U0);
params.pos.omega = fastopt.param('fix',v.pos.omega);
% Electrolyte.
params.const.psi = fastopt.param('fix',v.const.psi);
params.const.W = fastopt.param('fix',v.const.W);
params.pos.tauW = fastopt.param('fix',v.pos.tauW);
params.pos.kappa = fastopt.param('fix',v.pos.kappa);
params.eff.tauW = fastopt.param('fix',v.pos.tauW);
params.eff.kappa = fastopt.param('fix',v.pos.kappa);
% Porous electrode.
params.pos.Rdl = fastopt.param('fix',v.pos.Rdl);
params.pos.Cdl = fastopt.param('fix',v.pos.Cdl);
params.pos.nDL = fastopt.param('fix',v.pos.nDL);
params.pos.Rf = fastopt.param('fix',v.pos.Rf);
params.pos.k0Theta = fastopt.param('fix',v.pos.k0Theta);
params.pos.k0Linear = fastopt.param('fix',v.pos.k0Linear);
params.pos.DsTheta = fastopt.param('fix',v.pos.DsTheta);
params.pos.DsLinear = fastopt.param('fix',v.pos.DsLinear);
params.pos.nF = fastopt.param('fix',v.pos.nF);
params.pos.tauF = fastopt.param('fix',v.pos.tauF);
params.pos.sigma = fastopt.param('fix',v.pos.sigma);
% Lithium-metal electrode.
params.neg.Rdl = fastopt.param('fix',v.neg.Rdl);
params.neg.Cdl = fastopt.param('fix',v.neg.Cdl);
params.neg.nDL = fastopt.param('fix',v.neg.nDL);
params.neg.Rf = fastopt.param('fix',v.neg.Rf);
params.neg.k0 = fastopt.param('fix',v.neg.k0);

% Free paramters.
% Reaction symmetry factors.
params.neg.alpha = fastopt.param;
params.pos.alphaLinear = fastopt.param('len',length(v.pos.k0Theta));

modelspec = fastopt.modelspec(params, ...
    'tempsdegC',TdegC,'TrefdegC',labLinearFit.arg.TrefdegC);


% Define optimization bounds ----------------------------------------------
lb.neg.alpha = 0;
ub.neg.alpha = 1;
lb.pos.alphaLinear = zeros(length(v.pos.k0Theta),1); 
ub.pos.alphaLinear = ones(length(v.pos.k0Theta),1);

init.neg.alpha = 0.5;
init.pos.alphaLinear = 0.5*ones(length(v.pos.k0Theta),1);

% Perform regression ------------------------------------------------------

% Collect weight matrix.
weights = cell(multiplicity,1);
for m = 1:multiplicity
    weights{m} = ones(length(freqLab{m}),length(socPctTrue{m}));
    if ~isempty(arg.WeightFcn)
        for idxSOC = 1:length(socPctTrue{m})
            for idxFreq = 1:length(freqLab{m})
                weights{m}(idxFreq,idxSOC) = arg.WeightFcn( ...
                    freqLab{m}(idxFreq),socPctTrue{m}(idxSOC),TdegC(m));
            end % for
        end % for
    end % if
end % for

[plotInit, plotUpdate] = uiH2ImpedanceCallbacks( ...
    modelspec,Zlab,freqLab,socPctTrue,TdegC);
data = fastopt.uiparticleswarm(@cost,modelspec,init,lb,ub, ...
    'PlotInitializeFcn',plotInit, ...
    'PlotUpdateFcn',plotUpdate, ...
    'UseParallel',arg.UseParallel);

% Collect output data.
data.model = setCellParam(initialModel,data.values);
data.Zmodel = getH2Impedance(data.values,freqLab,socPctTrue,TdegC,ocpData);
data.Zlab = Zlab;
data.freq = freqLab;
data.socPctTrue = socPctTrue;
data.TdegC = TdegC;
data.arg = arg;
data.type__ = 'ParameterEstimate';
data.origin__ = 'fitNonlinearEIS';

function J = cost(model)
    modelVect = fastopt.splittemps(model,modelspec);

    % Ignore singular matrix warnings in solving NL model.
    warning('off','MATLAB:nearlySingularMatrix');

    J = 0;
    for km = 1:multiplicity
        % Calculate impedance predicted by the linear EIS model.
        Zmodel = getH2Impedance( ...
            modelVect(km),freqLab{km},socPctTrue{km},TdegC(km),ocpData{km});
        % Compute total residual between model impedance and measured
        % impedance across all spectra.
        J = J + sum(weights{km}.*(abs(Zmodel-Zlab{km})./abs(Zlab{km})).^2,'all');
    end

    % Re-enable singular matrix warning.
    warning('on','MATLAB:nearlySingularMatrix');
end % cost()

end % fitLinearEIS()