function spectra = loadLabNLEIS(expDir,varargin)
%LOADLABNLEIS Load (NL)EIS spectra collected in the laboratory by
%  a Gamry potentiostat.
%
% spectra = LOADLABNLEIS(expDir) loads the (NL)EIS spectra collected by
%   a Gramry potentiostat into a MATLAB data structure.
%   EXPDIR is the output directory produced  by the Gamry (NL)EIS sequence. 
%   SPECTRA is a structure containing the linear and second-harmonic
%   impedance spectra (see description below).
%
% spectra = LOADLABNLEIS(expDir,'MaxFreqBands',B) explicitly specifies the
%   maximum number of distinct frequency bands for both linear and nonlinear
%   EIS as B. The default is B=10.
%
% The output, SPECTRA, is a structure with the following fields:
%   lin.Z   : matrix of linear spectra (dim1=freq, dim2=SOC) [V/A]
%   lin.freq: cyclic frequency vector for linear spectra [Hz]
%   lin.ocv : matrix of mean OCV measurements (dim1=band, dim2=SOC) [V]
%   h2.Z    : matrix of second-harmonic spectra (dim1=freq, dim2=SOC) [V/A^2]
%   h2.freq : cyclic frequency vector for second-harmonic spectra [Hz]
%   h2.ocv  : matrix of mean OCV measurements (dim1=band, dim2=SOC) [V]
%   socPct  : SOC vector for both linear and second-harmonoc spectra [%]
%   soc0Pct : cell SOC at start of the experiment [%]
%   QdisAh  : vector specifying total amount of charge removed from cell
%             before EIS expirements at each SOC setpoint were run [Ah]
%   TdegC   : average temperature measured during discharge scripts [degC]
%   cellName: string identifying the cell
%   instrumentName: string identifying the Gamry potentiostat
%   linDriftCorrect: whether or not linear spectra is drift-corrected
%   experimentSpec: structure of experiment parameter values from 
%     the corresponding Parameters.mat file
%   param   : structure parameter values supplied to this function
%
% -- Changelog --
% 2023.05.15 | Created | Wesley Hileman <whileman@uccs.edu>

isinteger = @(x)floor(x)==x;
parser = inputParser;
parser.addRequired('expDir',@(x)ischar(x)||isstring(x));
parser.addParameter('MaxFreqBands',10,@(x)isscalar(x)&&isinteger(x)&&x>=1);
parser.parse(expDir,varargin{:});
p = parser.Results; % structure of validated arguments

expSpec = load(fullfile(expDir,'Parameters.mat'));

% Parse linear EIS --------------------------------------------------------

% Determine number of frequency bands included.
for j = 1:p.MaxFreqBands
    fieldname = sprintf('eis%d_filenames',j);
    if ~isfield(expSpec,fieldname)
        break;
    end
end % for
bandCount = j-1;

% Adjoin impedance data from each band.
clear bands;
for j = bandCount:-1:1  % iterate bands
    freq = double.empty(length(expSpec.soc),0);
    Z = double.empty(length(expSpec.soc),0);
    ocv = nan(length(expSpec.soc),1);
    fieldname = sprintf('eis%d_filenames',j);
    for k = 1:size(expSpec.discharge_filenames,1)  % iterate SOC setpoints
        filename = strtrim(expSpec.(fieldname)(k,:));
        filepath = fullfile(expDir,filename);
        if isfile(filepath)
            eis = loadGamryEIS(filepath);
            freq(k,1:length(eis.freq)) = eis.freq(:);
            Z(k,1:length(eis.freq)) = eis.Z(:);
            ocv(k) = mean(eis.tables.OCVCURVE.Vm);
        else
            Z(k,:) = NaN;
        end
    end
    bands(j).freq = freq(1,:); % assume same frequencies at each SOC!
    bands(j).Z = Z; 
    bands(j).ocv = ocv;
end
spectra.lin.Z = [bands.Z].';
spectra.lin.freq = [bands.freq].';
spectra.lin.ocv = [bands.ocv].';


% Parse nonlinear EIS -----------------------------------------------------

% Determine number of frequency bands included.
for j = 1:p.MaxFreqBands
    fieldname = sprintf('nleis%d_filenames',j);
    if ~isfield(expSpec,fieldname)
        break;
    end
end % for
bandCount = j-1;

% Adjoin impedance data from each band.
clear bands;
for j = bandCount:-1:1  % iterate bands
    freq = double.empty(length(expSpec.soc),0);
    Z = double.empty(length(expSpec.soc),0);
    ocv = nan(length(expSpec.soc),1);
    fieldname = sprintf('nleis%d_filenames',j);
    for k = 1:size(expSpec.discharge_filenames,1)  % iterate SOC setpoints
        filename = strtrim(expSpec.(fieldname)(k,:));
        filepath = fullfile(expDir,filename);
        if isfile(filepath)
            eis = loadGamryEIS(filepath);
            freq(k,1:length(eis.freq)) = eis.freq(:);
            Z(k,1:length(eis.freq)) = eis.Zhh(:,2);
            ocv(k) = mean(eis.tables.OCVCURVE.Vm);
        else
            Z(k,:) = NaN;
        end
    end
    bands(j).freq = freq(1,:); % assume same frequencies at each SOC!
    bands(j).Z = Z; 
    bands(j).ocv = ocv;
end
spectra.h2.Z = [bands.Z].';
spectra.h2.freq = [bands.freq].';
spectra.h2.ocv = [bands.ocv].';


% Parse discharge data ----------------------------------------------------

% Determine Ah discharged before running EIS at each SOC setpoint.
QdisAh = zeros(size(expSpec.soc));
meanTemp = nan(size(expSpec.soc));
for k = 1:size(expSpec.discharge_filenames,1)
    filename = strtrim(expSpec.discharge_filenames(k,:));
    filepath = fullfile(expDir,filename);
    if isfile(filepath)
        dis = loadGamryDTA(fullfile(expDir,filename));
        tab = dis.tables.CURVE;
        time = tab.T;
        iapp = -tab.Im; % Gamry uses opposite sign convention for iapp
        QdisAh(k) = trapz(time,iapp)/3600;
        meanTemp(k) = mean(tab.Temp(:));
    end % if
end % for
spectra.soc0Pct = expSpec.soc0;
spectra.socPct = expSpec.soc(:);
spectra.QdisAh = QdisAh(:);
spectra.TdegC = mean(meanTemp(~isnan(meanTemp)));
spectra.cellName = expSpec.cell;
spectra.instrumentName = expSpec.instrument;
spectra.linDriftCorrect = expSpec.enable_drift_correct;

spectra.experimentSpec = expSpec;
spectra.param = p;
spectra.origin__ = 'processLabNLEIS';

end


