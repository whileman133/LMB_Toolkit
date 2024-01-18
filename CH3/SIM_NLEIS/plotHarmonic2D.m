% plotHarmonicImpedance2D.m
%
% Plot linear and second-harmonic impedance of a 2D electrode (e.g., the
% lithium-metal electrode of an LMB cell).
%
% -- Changelog --
% 2024.01.10 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath(fullfile('..','..'));
TB.addpaths;

% Constants.
sensFreq = logspace(0.5,5,1000); % frequencies over which to evalulate Z [Hz]
simFreq = logspace(0.5,4,40);
sensBeta = (0.3:0.1:0.7)';      % values of beta to use in sens analysis
TdegC = 25;                     % temperature [degC]
I = 0.1:0.1:0.5;                % current for simulation study [A]

% Construct model of 2D electrode.
m2D.i0 = 0.2;    % exchange current [A]
m2D.Cdl = 5e-3;  % double-layer capacitance [F]
m2D.beta = 0.3;  % reaction symmetry factor [-]


% Sensitivity study -------------------------------------------------------
spec.defaults = m2D;
spec.singl.values.beta = sensBeta;
sensData = fastopt.runSensitivityStudy(spec,@(x)calcHarmonic2D(sensFreq,x,TdegC));
Z1 = [sensData.results.output.h1];
Z2 = [sensData.results.output.h2];

lab = arrayfun(@(x)sprintf('\\beta=%.2f',x),sensBeta,'UniformOutput',false);
col = cool(length(sensBeta));

figure; colororder(col);
plot(real(Z1),-imag(Z1));
xlabel('Re');
ylabel('-Im');
title('Z1');
legend(lab,'Location','best');
setAxesNyquist;
thesisFormat;
print(fullfile('plots','sens-Z1-beta'),'-depsc');
print(fullfile('plots','sens-Z1-beta'),'-dpng');

figure; colororder(col);
plot(real(Z2),-imag(Z2));
xlabel('Re');
ylabel('-Im');
title('Z2');
legend(lab,'Location','best');
setAxesNyquist;
thesisFormat;
print(fullfile('plots','sens-Z2-beta'),'-depsc');
print(fullfile('plots','sens-Z2-beta'),'-dpng');


% Simulation study --------------------------------------------------------

Zmodel = calcHarmonic2D(sensFreq,m2D,TdegC);
Zsim1 = zeros(length(simFreq),length(I));
Zsim2 = zeros(length(simFreq),length(I));
for k = 1:length(I)
    fprintf('Simulating @ I=%.3f\n',I(k));
    data = simHarmonic2D(simFreq,m2D,I(k),TdegC);
    Zsim1(:,k) = data.h1;
    Zsim2(:,k) = data.h2;
end
% Plot nonlinear time waveforms.
simHarmonic2D(1,m2D,0.5,TdegC,true);

lab = [{'Analytic Solution'} arrayfun(@(x)sprintf('ode23 $i_\\mathrm{app}=%.1f\\,\\mathrm{A}$',x),I,'UniformOutput',false)];
col = [0 0 0; cool(length(I))];

figure; colororder(col);
plot(real(Zmodel.h1),-imag(Zmodel.h1)); hold on;
plot(real(Zsim1),-imag(Zsim1),'d');
xlabel('Re');
ylabel('-Im');
title('Z1');
legend(lab,'Location','south','Interpreter','latex');
setAxesNyquist;
thesisFormat;
print(fullfile('plots','sim-Z1-I'),'-depsc');
print(fullfile('plots','sim-Z1-I'),'-dpng');

figure; colororder(col);
plot(real(Zmodel.h2),-imag(Zmodel.h2)); hold on;
plot(real(Zsim2),-imag(Zsim2),'d');
xlabel('Re');
ylabel('-Im');
title('Z2');
%legend(lab,'Location','best','Interpreter','latex');
setAxesNyquist;
thesisFormat;
print(fullfile('plots','sim-Z2-I'),'-depsc');
print(fullfile('plots','sim-Z2-I'),'-dpng');