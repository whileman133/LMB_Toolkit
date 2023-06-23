function result = linKK(f,z,varargin)
%LINKK Test impedance spectra for causality, time-invariance, and linearity.
%
% Validate an impedance spectrum against the Kramers–Kronig relations by
% fitting an LTI model comprising a number of series-connected Voigt 
% (parallel resistor-capacitor) to the spectrum.
% Implements an algorithm that automatically selectes the number of
% elements (M) to avoid under- and over-fitting.
%
% In circuit form, the model is represented as follows:
% 
%                                   C1                  CM
%                               |---||----|         |---||----|
% --\/\/\-----)()()(-----||-----|         |-- ... --|         |--
%    R0         L0       C0     |--\/\/\--|         |--\/\/\--|
%                    (optional)     R1                  RM
%                                          (M times)
% 
% where the elements may take on either positive or negative values
% (negative values do not have physical meaning, but still satisfy the
% the KK relations).
%
% result = LINKK(f,z) fits a model without series integrator to the
%  impedance spectrum specified by frequency vector F and complex impedance
%  vector Z. The number of elements M is automatically selected. RESULT is
%  a structure with the following fields:
%  
%    result.M          The number of series Voigt elements fit.
%    result.mu         Under/over-fitting factor for the final fit model.
%    result.Zmodel     Impedance predicted by fit model corresponding to
%                      frequencies in the vector F.
%    result.residuals  Complex residuals, normalized to the magnitude of
%                      the true measured impedance ||Z||.
%
% result = LINKK(...,'c',c) sets the threshold for the under/over-fitting 
%   factor (mu) to c instead of the default 0.8. c must be between 0 and 1
%   and determines the selection of the number of series RC pairs M.
%
% result = LINKK(...,'IncludeIntegrator',true) includes a series integrator
%   (capacitor) in the model for fitting the low-frequency response of a
%   battery cell. The default model does not include an integrator.
%
% results = LINKK(...,'maxM',maxM) constrains the maximum number of series
%   elements to MAXM instead of the default 100.
%
% Adapted to MATLAB from https://github.com/ECSHackWeek/impedance.py,
% which is based on the method presented in:
%   Schönleber, M. et al. A Method for Improving the Robustness of
%   linear Kramers-Kronig Validity Tests. 
%   Electrochimica Acta 131, 20–27 (2014)
%   doi: 10.1016/j.electacta.2014.01.034
%   https://doi.org/10.1016/j.electacta.2014.01.034
%
% -- Changelog --
% 2022.12.26 | Created | Wesley Hileman <whileman@uccs.edu>

p = inputParser;
p.addRequired('f',@(x)isvector(x)&&all(x>0));
p.addRequired('z',@(x)isvector(x));
p.addParameter('c',0.8,@isscalar);
p.addParameter('IncludeIntegrator',false,@(x)isscalar(x)&&islogical(x)&&0<=x&&x<=1);
p.addParameter('maxM',100,@(x)isscalar(x)&&isinteger(x)&&x>=1);
p.parse(f,z,varargin{:});
c = p.Results.c;
useInteg = p.Results.IncludeIntegrator;
maxM = p.Results.maxM;

tauL = 1/2/pi/max(f);
tauH = 1/2/pi/min(f);

for M = 1:maxM
    % Logarithmically space time-constants of the RC circuits.
    tau = logspace(log10(tauL),log10(tauH),M);
    result = fitKK_(f,z,tau,useInteg);
    if result.mu <= c
        break;
    end
end

result.M = M;

end

function result = fitKK_(f,z,tau,useInteg)
    tau = tau(:)';       % Force row vector.
    z = z(:); f = f(:);  % Force column vector.
    magZ = abs(z);
    omega = 2*pi*f;
    s = 1j*omega;
    M = length(tau);
    coeffRn = 1./(1+s*tau);

    % The assumed parameter vector is: x = [R0, L0, 1/C0, R1, R2, ... RM]'.

    % Build measurement vector (y).
    y = [real(z)./magZ; imag(z)./magZ];

    % Build measurement matrix (H). Columns correspond to elements
    % in the circuit model.
    if useInteg
        Hrealz = zeros(length(z),3+M);   Himagz = zeros(length(z),3+M);
        Hrealz(:,1) = 1;                 Himagz(:,1) = 0;           % R0
        Hrealz(:,2) = 0;                 Himagz(:,2) = omega;       % L0
        Hrealz(:,3) = 0;                 Himagz(:,3) = -1./omega;   % 1/C0
        Hrealz(:,4:end) = real(coeffRn); Himagz(:,4:end) = imag(coeffRn); % R1...RM
        H = [Hrealz./magZ; Himagz./magZ];
    else
        Hrealz = zeros(length(z),2+M);   Himagz = zeros(length(z),2+M);
        Hrealz(:,1) = 1;                 Himagz(:,1) = 0;           % R0
        Hrealz(:,2) = 0;                 Himagz(:,2) = omega;       % L0
        Hrealz(:,3:end) = real(coeffRn); Himagz(:,3:end) = imag(coeffRn); % R1...RM
        H = [Hrealz./magZ; Himagz./magZ];
    end

    % Find parameter vector by least-squares regression.
    x = H\y;

    % Compute impedance predicted by the fit model.
    if useInteg
        elem.R0 = x(1);
        elem.L0 = x(2);
        elem.C0inv = x(3);
        elem.Rn = x(4:end);
        Zmodel = elem.R0 + s*elem.L0 + elem.C0inv./s + coeffRn*elem.Rn;
    else
        elem.R0 = x(1);
        elem.L0 = x(2);
        elem.Rn = x(3:end);
        Zmodel = elem.R0 + s*elem.L0 + coeffRn*elem.Rn;
    end

    % Compute under/over-fitting factor and complex residuals.
    Rn = elem.Rn;
    mu = 1 - sum(abs(Rn(Rn<0)))/sum(abs(Rn(Rn>=0)));
    if isnan(mu) || mu < 0
        mu = 1;
    end
    residuals = (Zmodel-z)./abs(z);

    % Store result in output structure.
    result.mu = mu;
    result.Zmodel = Zmodel;
    result.residuals = residuals;
    result.elements = elem;
end