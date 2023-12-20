function hcycData = processSimHalfCycle(simData)
%PROCESSSIMHALFCYCLE Extract OCP and resistance from half-cycle discharge
%  data from simHalfCycle.m

% Fetch time, iapp(t), vcell(t) vectors.
halfcycData = [simData.series];
Z0 = [simData.series.z0].';  % initial SOC [fractional]
Iapp = [halfcycData.iapp];
Vcell = [halfcycData.vcell];
time = zeros(size(Iapp));
for k = 1:size(time,2)
    time(:,k) = halfcycData(k).time(:);
end

% Compute average SOC of porous electrode as a function of time.
Q = simData.Q;
Zavg = zeros(size(time)); % average SOC vs time [fractional]
for k = 1:size(Zavg,2)
    Zavg(:,k) = Z0(k) - cumtrapz(time(:,k),Iapp(:,k))/Q/3600;
end

% Get Vcell, Iapp over common SOC.
zavg = linspace(min(Zavg,[],'all'),max(Zavg,[],'all'),1000).';
IappNorm = zeros(length(zavg),size(Iapp,2));
VcellNorm = zeros(length(zavg),size(Vcell,2));
for k = 1:size(Iapp,2)
    IappNorm(:,k) = interp1(Zavg(:,k),Iapp(:,k),zavg,'linear','extrap');
    VcellNorm(:,k) = interp1(Zavg(:,k),Vcell(:,k),zavg,'linear','extrap');
end

% Estimate OCP and resistance at each lithiation point using least-squares
% linear regression.
UocpEst = zeros(size(zavg));
RdcEst = zeros(size(zavg));
for k = 1:size(IappNorm,1)
    % Collect all current magnitudes applied and cell voltages measured 
    % at this SOC point into column vectors.
    iapp = IappNorm(k,:).';
    vcell = VcellNorm(k,:).';
    % Compute least-squares solution to the linear system: 
    % vcell(z) = Uocp(z) - Rcell(z)*iapp(z)
    Ieqls0 = abs(iapp)<Q/200;  % logical indicies to zero current
    if sum(~Ieqls0)<2
        % Not enough nonzero current data-points to compute Uocp and Rcell
        % using LS; approximate Uocp only.
        UocpEst(k) = mean(vcell(Ieqls0));
        RdcEst(k) = NaN;
        continue;
    end
    H = [ones(size(iapp)) -iapp];  % measurement matrix
    x = H\vcell;
    UocpEst(k) = x(1);
    RdcEst(k) = x(2);
end

% Collect output data.
hcycData.ocp.socAvgPct = zavg*100;
hcycData.ocp.Uocp = UocpEst;
hcycData.ocp.TdegC = simData.TdegC;
% Remove NaNs from dc resistance estimate!
hcycData.res.socAvgPct = zavg*100;
hcycData.res.Rdc = RdcEst;
hcycData.res.TdegC = simData.TdegC;
hcycData.VcellNorm = VcellNorm;
hcycData.IappNorm = IappNorm;
hcycData.simData = simData;

end