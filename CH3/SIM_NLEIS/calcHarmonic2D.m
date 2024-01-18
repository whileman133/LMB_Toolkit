function data = calcHarmonic2D(freq,model,TdegC)
%CALCHARMONIC2D Calculate linear and 2nd harmonic impedance of 2D electrode.
%
% -- Changelog --
% 2024.01.10 | Created | Wesley Hileman <whileman@uccs.edu>

omega = 2*pi*freq(:);
f = TB.const.f(TdegC);
F = TB.const.F;
R = TB.const.R;
T = TdegC+273.15;
Rct = 1/model.i0/f;
Cdl = model.Cdl;
b = model.beta;

data.h1 = Rct./(1+1j*omega*Rct*Cdl);
data.h2 = -((1-b)^2-b^2)*F*Rct*Rct/4/R/T...
    ./(1+1j*2*omega*Rct*Cdl)...
    ./(1+1j*1*omega*Rct*Cdl).^2;

end