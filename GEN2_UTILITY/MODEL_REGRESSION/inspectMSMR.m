function [U0, X, omega] = inspectMSMR(ocp, J, varargin)
%INSPECTMSMR Roughly estimate MSMR parameters from lab data over
%  relative composition (electode SOC). Plots differential capacity
%  for manual inspection.
%
% This works well only when all of the galleries are observable
% in the OCP curve. Manual fudging of this guess is required 
% for materials such as NMC.
%
% The reported X values should be scaled by (theta_max-theta_min).
% Manually adjusting the X is neccesary to ensure sum(X)=1 
% after scaling. In particular, consider stretching the X of 
% the galleries at the "edges" of the differential capacity 
% curves whose X may be larger than shown in collected data.
%
% -- Usage --
% [U0,X,omega] = inspectMSMR(ocp,J)
%
% -- Input --
%  ocp   = struct or struct array of OCP data (see fields below)
%  J     = number of galleries (actual found may be fewer)
%
% -- Output --
%   U0, X, omega = vectors of rough estimates of MSMR parameters
%
% -- Changelog --
% 2024.01.18 | Created | Wesley Hileman <whileman@uccs.edu>

% Constants.
R = 8.3144598;      % Molar gas constant [J/mol K]
F = 96485.3329;     % Faraday constant [C/mol]

% Parse arguments.
parser = inputParser;
parser.addRequired('ocp',@(x)isscalar(x)&&isstruct(x));
parser.addRequired('J',@(x)isscalar(x)&&x>0);
parser.parse(ocp,J,varargin{:});
arg = parser.Results;  % struct of validated arguments

% Collect arguments.
J = arg.J;

% Collect voltage-v-SOC curve.
TdegC = ocp.TdegC;
Z = ocp.Z;
V = ocp.U;
dZ = abs(ocp.dZ);
d2Z = calc_dZ(dZ,V);
d3Z = calc_dZ(d2Z,V);

% Normalize.
d2Z = d2Z/max(abs(d2Z));
d3Z = d3Z/max(abs(d3Z));

% Locate stationary points of dZ within an approximate margin.
[d2ZVal,ind_d2Zeqls0] = findpeaks(-abs(d2Z),'SortStr','descend','NPeaks',2*J);
ind_d2Zeqls0 = ind_d2Zeqls0(d2ZVal>-0.2);
% Separate maxima from minima.
d3ZVal = d3Z(ind_d2Zeqls0);
ind_d2Zmax = ind_d2Zeqls0(d3ZVal<+0.01);  % maxima - concave down
ind_d2Zmin = ind_d2Zeqls0(d3ZVal>-0.01);  % minima - concave up
% Limit to J maxima, J-1 minima at most.
limmax = min(J,length(ind_d2Zmax));
limmin = min(J-1,length(ind_d2Zmin));
ind_d2Zmax = ind_d2Zmax(1:limmax);
ind_d2Zmin = ind_d2Zmin(1:limmin);

% Estimate U0 and X
[U0, indsort] = sort(V(ind_d2Zmax),'descend');
dZU0 = dZ(ind_d2Zmax);
dZU0 = dZU0(indsort);
U0 = U0(:);
dZU0 = dZU0(:);
dX = sort(Z(ind_d2Zmin));
X  = diff([0; dX(:); 1]);

% Check that we got the right number of U0 and X values.
if length(U0) < J
    fprintf(['WARNING: number of U0 found (%d) is less than ' ...
        'target number of phases (%d).\n'],length(U0),J);
end
if length(X) < J
    fprintf(['WARNING: number of X found (%d) is less than ' ...
        'target number of phases (%d).\n'],length(X),J);
end

% Estimate omega if U0 and X are consistent.
if length(X) ~= length(U0)
    omega = []; 
    fprintf(['WARNING: X-U0 mismatch. Could not estimate omega. ' ...
        'Incomplete MSMR estimate. ']);
else
    f = F/R/(TdegC+273.15);
    omega = f*X./dZU0/4;
end

% Plotting ----------------------------------------------------------------
if ~isempty(omega)
    % Compute and plot MSMR estimate.
    params.U0 = U0;
    params.X = X;
    params.omega = omega;
    params.thetamin = 0;
    params.thetamax = 1;
    ocpData = MSMR(params).ocp('voltage',V,'TdegC',TdegC,'npoints',100);
end

% OCP-vs-SOC
figure;
plot(Z,V); hold on;
if ~isempty(omega)
    plot(ocpData.theta,ocpData.Uocp,':');
end
yline(U0,'k');
xlabel('SOC');
ylabel('Uocp');
title('inspectMSMR: Uocp-vs-SOC curve');
if ~isempty(omega)
    legend('Lab','Guess');
end
thesisFormat;

% Differential-Capaity-vs-OCP
figure;
plot(V,dZ); hold on;
if ~isempty(omega)
    plot(ocpData.Uocp,abs(1./ocpData.dUocp),':');
end
xline(U0);
xlabel('Uocp');
ylabel('dZ/dU');
title('inspectMSMR: dZ-vs-Uocp curve');
if ~isempty(omega)
    legend('Lab','Guess');
end
thesisFormat;

% Differential-Capacity-vs-SOC
figure;
plot(Z,dZ); hold on;
if ~isempty(omega)
    plot(ocpData.theta,abs(1./ocpData.dUocp),':');
end
xlabel('SOC');
ylabel('dZ/dU');
title('inspectMSMR: dZ-vs-SOC curve');
if ~isempty(omega)
    legend('Lab','Guess');
end
thesisFormat;

end

function dZ = calc_dZ(Z,V)
  dZ = diff(Z(:))./diff(V(:));
  dZ = ([dZ(1);dZ]+[dZ;dZ(end)])/2; % Avg of fwd/bkwd diffs
end