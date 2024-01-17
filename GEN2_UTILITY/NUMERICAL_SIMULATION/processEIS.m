function spectra = processEIS(simData,varargin)
%PROCESSEIS Process time-domain data collected by EIS simulation in COMSOL.
%
% spectra = PROCESSEIS(simData) converts the time-domain data collected by
%   COMSOL EIS simulation performed using SIMFOM.m into frequency-domain
%   "transfer-function" quantities. SIMDATA is the output produced by
%   SIMFOM.m. SPECTRA is a structure containing the transfer-function
%   spectra (see description below), which is obtained using 
%   the fast-Fourier-transform (FFT) algorithm.
%
% spectra = PROCESSEIS(...,'NumHarmonics',H) specifies the number of
%   harmonic spectra to compute. H=1 computes the linear spectra only;
%   H=2 computes both linear and second-harmonic spectra; H=3 computes
%   linear, second-, and third-harmonic spectra. Default H=2.
%
% spectra = PROCESSEIS(...,'EvalLinTF',true) also evalulates the linear
%   transfer-function model for the cell and stores the result in the
%   output for comparison to the linear COMSOL spectra.
%
% spectra = PROCESSEIS(...,'EvalLinTF',true,'NumTFFreqPoints',N) also 
%   specifies the number of frequency-points at which to evalulate the 
%   linear transfer-functions. Default N=1000.
%
% The output, SPECTRA, is a structure with the following fields:
%   lin   : structure array of linear spectra for each variable [VariableUnit/A]
%   h2    : (if H>=2) structure array of second-harmonic spectra [VariableUnit/A^2]
%   h3    : (if H>=3) structure array of third-harmonic spectra [VariableUnit/A^3]
%   freq  : frequency vector corresponding to entries in the structure arrays [Hz]
%   xlocs : structure mapping electrochemical variables to the x-locations
%           where the spectra are evalulated (only for internal cell variables)
%   tf    : (present if EvalLinTF==true) structure array of linear spectra 
%           evalulated from transfer functions for each variable.
%   tfFreq: (present if EvalLinTF==true) frequency vector corresponding to 
%           the tf structure array [Hz].
%   cellModel: cell models structure from which impedance was simulated
%   TdegC : temperature at which simulation was run [degC].
%   socPct: vector of SOC setpoints at which simulation was run [%].
%   param : structure of parameter values supplied to the function.
% 
% -- Examples --
% SINGLE SOC SETPOINT
% Let:
%   spectra = processEIS(simData);
% where simData.arg.socPct was set to 50%.
% 
% 1. Fetch the linear impedance spectrum:
%    Z1 = [spectra.lin.Zcell];  % row vector, dim1=freq
%
% 2. Fetch the linear Thetae(jw)/Iapp(jw) spectrum at x=3:
%    ind = find(spectra.xlocs.Thetae==3,1,'first')
%    Thetae = [spectra.lin.Thetae]; % matrix, dim1=x-loc, dim2=freq
%    Thetae3 = Thetae(ind,:); % row vector, dim1=freq
%
% MULTIPLE SOC SETPOINTS
% Let:
%   spectra = processEIS(simData);
% where simData.arg.socPct was set to [5 95]%.
%
% 1. Fetch linear impedance spectra at all SOC setpoints:
%    Z1 = reshape([spectra.lin.Zcell],size(spectra.lin)]; % matrix, dim1=soc, dim2=freq
%
% 2. Fetch the linear Thetae(jw)/Iapp(jw) spectrum at x=3 and all SOCs:
%    ind = find(spectra.xlocs.Thetae==3,1,'first')
%    Thetae = [spectra.lin.Thetae]; % matrix, dim1=x-loc, dim2=freq(repeat for each soc)
%    Thetae3 = reshape(Thetae(ind,:),size(spectra.lin)); % matrix, dim1=soc, dim2=freq
%
% -- Changelog --
% 2023.12.27 | Add Eta TF | Wesley Hileman
% 2023.06.11 | Support multiple SOC setpoints | Wesley Hileman
% 2023.04.05 | Created | Wesley Hileman <whileman@uccs.edu>

isinteger = @(x)floor(x)==x;
parser = inputParser;
parser.addRequired('simData',@(x)isstruct(x)&&strcmp(x.origin__,'simEIS'));
parser.addParameter('NumHarmonics',2,@(x)isscalar(x)&&isinteger(x)&&x>=1);
parser.addParameter('EvalLinTF',false,@(x)isscalar(x)&&islogical(x));
parser.addParameter('NumTFFreqPoints',1000,@(x)isscalar(x)&&isinteger(x));
parser.parse(simData,varargin{:});
p = parser.Results; % structure of validated arguments
H = p.NumHarmonics; % number of harmonics to retain

varNames = setdiff(fieldnames(simData.ss),{'param','time'});
harmNames = arrayfun(@(x)sprintf('h%d',x),1:H,'UniformOutput',false);
harmNames{1} = 'lin';

for h = 1:H
    spectra.(harmNames{h}) = [];
end

for kz = size(simData.ss,1):-1:1
    for kf = size(simData.ss,2):-1:1
        freqkData = simData.ss(kz,kf);
        socPct = simData.arg.socPct(kz);
    
        % Construct FFT frequency vector.
        fs = freqkData.param.fs;  % sampling rate [Sa/s]
        N = freqkData.param.N;    % number of samples available for fft
        f0 = freqkData.param.f0;  % fundemental frequency [Hz]
        tfFreq = 0:(fs/N):(fs/2); % frequency vector [Hz]
        fh = (1:H)*f0;            % harmonic frequencies [Hz]       
        [resid,ih] = min(abs(fh-tfFreq')); % harmonic indicies into fft vector
        if any(resid./fh>1e-6)
            warning('Possible FFT leakage for SOC=%.3f%% f=%.5g!',socPct,f0);
        end
    
        % Compute frequency spectra of the electrochemical variables.
        Vars = struct;
        for j = 1:length(varNames)
            varName = varNames{j};
            timePosData = freqkData.(varName);
            freqPosData = fft(timePosData);
            % Use single-sided spectrum.
            freqPosData = 2*freqPosData(1:length(tfFreq),:)/size(freqPosData,1);
            % Retain only the 1:H harmonics.
            Vars.(varName) = freqPosData(ih,:);
        end
    
        % Normalize to input current phasor to obtain "transfer function" form.
        II = Vars.Iapp(1).^(1:H).';
        for j = 1:length(varNames)
            varName = varNames{j};
            Vars.(varName) = Vars.(varName)./II;
        end
    
        % Compute cell impedance.
        Vars.Zcell = -Vars.Vcell;  % Don't forget the negative sign!
    
        % Store results.
        varsn = fieldnames(Vars);
        for j = 1:length(varsn)
            varName = varsn{j};
            freqPosData = Vars.(varName);
            for h = 1:H
                % Use column vector for variables with more than 
                % one x-location so that it's easier to collect results using
                % the structure-array notation [structureArray.fieldName].
                spectra.(harmNames{h})(kz,kf).(varName) = freqPosData(h,:).';
            end
        end
    end % for freq
end % for soc

spectra.freq  = simData.arg.freq;
if isfield(simData,'xlocs')
    spectra.xlocs = simData.xlocs;
else
    % An earlier version of simEIS did not include the xlocs field.
    spectra.xlocs = simData.arg.Vars;
end
spectra.xlocs.Iapp  = [];
spectra.xlocs.Vcell = [];
spectra.xlocs.Zcell = [];

if p.EvalLinTF
    % Compare linear COMSOL spectra to those predicted by transfer functions.

    % Structure mapping COMSOL variables to TF functions.
    % See end of this file for definitions of custom functions!
    var2tf.Thetae = @tfThetae;
    var2tf.Idl = @tfIdl;
    var2tf.If = @tfIf;
    var2tf.Ifdl = @tfIfdl;
    var2tf.Phie = @getPhieTF;    % ground ref. at x=0- (at neg cc)
    var2tf.PhieTilde = @tfPhie;  % ground ref. at x=0+ (in electrolyte)
    var2tf.PhisTilde = @tfPhis;
    var2tf.Phise = @tfPhiseInt;
    var2tf.Thetass = @tfThetassInt;
    var2tf.Eta = @tfEta;
    var2tf.Vcell = @getVcellTF;
    var2tf.Zcell = @getZcellTF;
    
    cellModel = convertCellModel(simData.arg.cellModel,'LLPM');
    socPct = simData.arg.socPct;
    TdegC = simData.arg.TdegC;
    tfFreq = logspace( ...
        log10(min(spectra.freq)),log10(max(spectra.freq)),p.NumTFFreqPoints);
    ss = 1j*2*pi*tfFreq;
    
    varNames = setdiff(fieldnames(spectra.lin),{'Iapp'});
    spectra.tf = [];
    for kz = length(socPct):-1:1
        mod = evalSetpoint(cellModel,ss,socPct(kz)/100,TdegC+273.15);
        for kv = 1:length(varNames)
            varName = varNames{kv};
            if ~isfield(var2tf,varName)
                continue;
            end
            xlocs = spectra.xlocs.(varName);
            tfValues = var2tf.(varName)(ss,xlocs,mod);
            for kf = length(ss):-1:1
                spectra.tf(kf,kz).(varName) = tfValues(:,kf);
            end % for freq
        end % for var
    end % for soc
    
    spectra.tfFreq = tfFreq;
end

spectra.cellModel = simData.arg.cellModel;
spectra.TdegC = simData.arg.TdegC;
spectra.socPct = simData.arg.socPct;
spectra.param = p;
spectra.origin__ = 'processEIS';

end

function PhieTF = getPhieTF(s,x,cellModel)
    % Phie with ground reference at x=0- (at negative electrode 
    % current-collector).
    Phise0 = tfPhiseInt(s,0,cellModel);
    PhieTilde = tfPhie(s,x,cellModel);
    PhieTF = PhieTilde - Phise0;
end

function VcellTF = getVcellTF(s,~,cellModel)
    Phise = tfPhiseInt(s,[0 3],cellModel);
    PhieTilde3 = tfPhie(s,3,cellModel);
    Rc = getCellParams(cellModel,'const.Rc');
    VcellTF = -Phise(1,:) + Phise(2,:) + PhieTilde3 - Rc;
end

function ZcellTF = getZcellTF(s,~,cellModel)
    % Negate cell-voltage TF!
    ZcellTF = -getVcellTF(s,[],cellModel);
end
