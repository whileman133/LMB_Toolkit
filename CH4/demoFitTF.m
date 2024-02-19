% demoFitTF.m
%
% Verify the TF model regression using synthetic data.
%
% -- Changelog --
% 2024.02.03 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath('..');
TB.addpaths;
rng(0);  % make results repeatable

% Constants.
synthCell = 'cellLMO-P2DM.xlsx';
useParallel = true; % run regression using parallel processing
deltaPct = 50; % amount by which initial estimates can differ from truth [%] 
TdegC = [15 25 40]; % temperature setpoints [degC]
socPct = 100:-5:5; % SOC setpoints [%]
freq = logspace(-4,5,100); % frequency points [Hz]
s = 1j*2*pi*freq;

% Load model of true cell.
p2dm = loadCellModel(synthCell);
rlwrm = convertCellModel(p2dm,'RLWRM');
% Use higher W so that electrolyte contributes significantly to impedance.
ptmp.const.W = 2;
rlwrm = setCellParam(rlwrm,ptmp);

% Setup fixed parameters.
truth = getCellParams(rlwrm);
fix.neg.Rf  = 0;                  % Lump with package resistance
fix.neg.Rdl = 0;                  % Rdl<<Rct, so assume Rdl~=0 for neg
fix.pos.sigma = truth.pos.sigma;  % TF model very insensitive to sigma
% Since nDL=nF=1, corresponding time constants are meaningless.
fix.neg.tauDL = truth.neg.tauDL;
fix.pos.tauDL = truth.pos.tauDL;
fix.pos.tauF = truth.pos.tauF;

% Fetch true OCP parameters (assume already known for TF regression).
ocp = getCellParams( ...
    rlwrm,'const.Q pos.U0 pos.X pos.omega pos.theta0 pos.theta100');
ocp.pos.Q = ocp.const.Q;
ocp = ocp.pos;

% Simulate EIS measurements by evalulating TF model.
clear eis;
pkg = getCellParams(rlwrm,'pkg.*');
for kt = length(TdegC):-1:1
    Z = zeros(length(freq),length(socPct));
    for kz = length(socPct):-1:1
        tfData = tfLMB( ...
            s,rlwrm,'TdegC',TdegC(kt),'socPct',socPct(kz));
        Z(:,kz) = tfData.h11.tfZcell();
        % !!! Important: add contribution of cell package!
        Z(:,kz) = Z(:,kz) + pkg.R0 + pkg.L0*s(:);
    end
    tmp.socPct = socPct(:);
    tmp.freq = freq(:);
    tmp.Z = Z;
    tmp.TdegC = TdegC(kt);
    eis(kt) = tmp;
end

% Generate initial guess - skew the truth a bit.
exclude = {...
    'alpha',...
    'U0','X','omega','theta0','theta100',...
    'vmin','vmax','Tref',...
    'brugDeKappa','brugSigma',...
};  % names of parameter values not to randomize
p = getCellParams(p2dm);
regs = fieldnames(p);
for kr = 1:length(regs)
    reg = regs{kr};
    pnames = fieldnames(p.(reg));
    for kp = 1:length(pnames)
        pname = pnames{kp};
        if any(strcmpi(pname,exclude))
            continue;
        end
        val = p.(reg).(pname);
        eps = (deltaPct/100)*(2*rand(size(val))-1);
        val = val.*(1+eps);
        if any(strcmpi(pname,{'nDL','nF','sEps','eEps'}))
            % Ensure CPE exponents / volume fractions remain between 0 and 1!
            val(val<0) = 0;
            val(val>1) = 1;
        end
        p.(reg).(pname) = val;
    end
end
init = setCellParam(p2dm,p);

% Regress TF model.
regspec = regTF('msmr',eis,ocp,init,'fix',fix);
fitdata = fitTF(regspec,'UseParallel',useParallel);

% Save results.
save('demoFitTF.mat','init','regspec','fitdata','p2dm','rlwrm');