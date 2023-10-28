function Z = getLinearImpedance(model,freq,socPct,TdegC,ocpData)
%GETLINEARIMPEDANCE Calculate TF impedance at given frequency / SOC points.
%
% -- Usage --
% Z = getLinearImpedance(cellParams,freq,socPct,TdegC) calculates the
%   linear impedance of a cell at the cyclic frequencies FREQ and soc
%   setpoints SOCPCT given the structure of parameter values CELLPARAMS.
%   The computation is performed at temperature TDEGC.
%
% Z = getLinearImpedance(cellModel,freq,socPct,TdegC) performs the same
%   calculation given the full cell model instead of the cell parameters.
%
% Z = getLinearImpedance(...,ocpData) uses the cached OCP data OCPDATA,
%   corresponds to the soc setpoints in SOCPCT, to speed up the impedance
%   computation.

nfreq = length(freq);
nsoc = length(socPct);
s = 1j*2*pi*freq;  % Laplace variable

if isCellModel(model)
    params = getCellParams(model,'TdegC',TdegC);
else
    params = model;
end

useCachedOCP = false;
if exist('ocpData','var') && ~isempty(ocpData)
    ocpModel = MSMR(ocpData);
    ctData = ocpModel.RctCachedOCP(params.pos,ocpData);
    dsData = ocpModel.DsCachedOCP(params.pos,ocpData);
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