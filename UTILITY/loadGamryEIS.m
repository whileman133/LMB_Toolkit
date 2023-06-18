function eis = loadGamryEIS(filename,varargin)
%LOADGAMRYEIS Extract EIS information from a Gamry .DTA file.

% Number of harmonics appearing in the DTA file.
nharm = 10;

eis = loadGamryDTA(filename,varargin{:});
ztab = eis.tables.ZCURVE;
eis.freq = ztab.Freq;

% Check for drift-corrected data.
if any(strcmpi("ZrealDrCor",ztab.Properties.VariableNames))
    eis.Z = ztab.ZrealDrCor+1j*ztab.ZimagDrCor;
else
    eis.Z = ztab.Zreal+1j*ztab.Zimag;
end

% Check for harmonic-response data.
if any(strcmpi("Vthd",ztab.Properties.VariableNames))
    eis.Vthd = ztab.Vthd;
    eis.Ithd = ztab.Ithd;
    Vh = zeros(length(eis.freq),nharm);
    Ih = Vh;
    for h = 1:nharm
        nameVR = sprintf('VHr%d',h); nameIR = sprintf('IHr%d',h);
        nameVI = sprintf('VHi%d',h); nameII = sprintf('IHi%d',h);
        Vh(:,h) = ztab.(nameVR)+1j*ztab.(nameVI);
        Ih(:,h) = ztab.(nameIR)+1j*ztab.(nameII);
    end
    eis.Vh = Vh;
    eis.Ih = Ih;
    % Note: Gamry current convention is (+) for charge current,
    % so we do not need a negative sign here.
    eis.Zhh = Vh./(Ih(:,1).^(1:nharm));
end

end