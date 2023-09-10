function Z = getH2Impedance(model,freq,socPct,TdegC,ocpData)
%GETH2IMPEDANCE Calculate 2nd-harmonic impedance at given frequency / 
%   SOC points.
%
% -- Usage --
% Z = getH2Impedance(cellParams,freq,socPct,TdegC) calculates the
%   second-harmonic impedance of a cell at the cyclic frequencies FREQ 
%   and soc setpoints SOCPCT given the structure of parameter values 
%   CELLPARAMS. The computation is performed at temperature TDEGC.
%
% Z = getH2Impedance(cellModel,freq,socPct,TdegC) performs the same
%   calculation given the full cell model instead of the cell parameters.
%
% Z = getH2Impedance(...,ocpData) uses the cached OCP data OCPDATA,
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
if exist('ocpData','var')
    % Compute Rct and Ds from cached OCP data.
    ocpModel = MSMR(ocpData);
    ctData = ocpModel.RctCachedOCP(params.pos,ocpData);
    dsData = ocpModel.DsCachedOCP(params.pos,ocpData);
    useCachedOCP = true;
end

% Calculate impedance at each SOC setpoint.
Z = zeros(nfreq,nsoc);
for k = 1:nsoc
    if useCachedOCP
        % Update OCP, charge-transfer resistance, and solid diffusivity 
        % at each SOC setpoint manually for improved performance.
        params.pos.Uocp = ocpData.Uocp(k);
        params.pos.dUocp = ocpData.dUocp(k);
        params.pos.d2Uocp = ocpData.d2Uocp(k);
        params.pos.Rct = ctData.Rct(k);
        params.pos.Rct2Inv = ctData.Rct2Inv(k);
        params.pos.Ds = dsData.Ds(k);
    end
    tfData = tfLMB(s,params,'TdegC',TdegC,'socPct',socPct(k),'Calc22',true);
    Z(:,k) = tfData.h22.tfVcell();
end % for
end