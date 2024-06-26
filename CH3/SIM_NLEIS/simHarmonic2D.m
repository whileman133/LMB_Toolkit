function data = simHarmonic2D(freqHz,model,I,TdegC,makePlots)
%SIMHARMONIC2D Determine linear and 2nd and 3rd order nonlinear
%  impedance spectra of a two-dimensional electrode by simulation.
%
% The electrode-electrolyte interface is modeled by the Butler-Volmer 
% equation for reaction flux in parallel with a simple double-layer 
% capacitance.
%
% There are three model parameters:
%   i0 [A]    : The exchange current, assumed constant.
%   beta [-]  : The symmetry factor of the reduction reaction.
%   Cdl [F]   : The double-layer capacitance.
%
% data = sim2DElectrodeNLEIS(freqHz,model,TdegC) determines the 1st through 
%   3rd order impedance spectra of the electrode|electrolyte junction at
%   parameter values specified by the strutcure MODEL at frequencies 
%   FREQHZ. DATA is a structure with fields h1, h2, and h3 corresponding
%   to the 1st through 3rd harmonic impedance, respectively. 
%   The temperature in degC is set by TdegC. If any model paramter values
%   are vectors, a corresponding number of columns is returned for each
%   parameter value for each harmonic in the impedance structure Z.

if ~exist('makePlots','var')
    makePlots = false;
end

F = 9.6485e+04;     % Faraday constant [C/mol]
R = 8.31446262;     % Molar gas constant [J/mol·K]
T = TdegC+273.15;   % Absolute temperature [K]

i0 = model.i0(:)';      % Exchange current [A]
beta = model.beta(:)';  % Symmetry factor of the reduction reaction
Cdl = model.Cdl(:)';    % Double-layer capacitance [F]
nparamsets = max([length(i0) length(beta) length(Cdl)]);
if length(i0)==1
    i0 = i0*ones(1,nparamsets);
end
if length(beta)==1
    beta = beta*ones(1,nparamsets);
end
if length(Cdl)==1
    Cdl = Cdl*ones(1,nparamsets);
end

ff = freqHz(:);  % Force column vector!
f = F/R/T;

phise1 = zeros(length(ff),nparamsets);
phise2 = phise1;
phise3 = phise1;
for kp = 1:nparamsets
    for kf = 1:length(ff)
        fprintf('.');
        data = sim(ff(kf),i0(kp),beta(kp),Cdl(kp));
        phise1(kf,kp) = data.Phiseh(1);
        phise2(kf,kp) = data.Phiseh(2);
        phise3(kf,kp) = data.Phiseh(3);
        if kf == 1 && makePlots
            % Yes, this is bad practice, but in a crunch it's fast and
            % easy.
            figure;
            plot(data.time,data.iapp);
            title('iapp');
            xlabel('time');
            ylabel('iappt');
            thesisFormat();
            print(fullfile('plots','sim-time-iapp'),'-depsc');
            print(fullfile('plots','sim-time-iapp'),'-dpng');
            figure;
            plot(data.time,data.phise);
            title('phise');
            xlabel('time');
            ylabel('phiset');
            thesisFormat();
            print(fullfile('plots','sim-time-phise'),'-depsc');
            print(fullfile('plots','sim-time-phise'),'-dpng');
            figure;
            stem(data.freq,abs(data.Phise),'filled','MarkerSize',10);
            title('MagSpectrum');
            xlabel('Freq');
            ylabel('Mag');
            thesisFormat();
            print(fullfile('plots','sim-fft'),'-depsc');
            print(fullfile('plots','sim-fft'),'-dpng');
        end
    end
    fprintf('\n');
end

data.h1 = phise1;
data.h2 = phise2;
data.h3 = phise3;


function data = sim(f0,i0,beta,Cdl)
    Rct = 1./i0/f;
    fs = f0*100;                    % Sampling rate [Sa/s].
    t0 = 10*Rct*Cdl;                % Transient interval [s].
    tss = 10/f0;                    % Steady-state interval [s].
    t0 = t0+1/f0-rem(t0,1/f0);      % (round transient interval to start of next cycle of input)
    tt = 0:(1/fs):(t0+tss);         % Time vector [s].
    ss = tt>t0; % logical indicies to steady-state interval
    iapp = I*sin(2*pi*f0*tt);
    [time,phise] = ode23(@(t,y)getRHS(t,y,f0,i0,beta,Cdl),tt,0);

    % Compute frequency spectrum.
    Phise = fft(phise(ss));
    freq = 0:(fs/length(Phise)):(fs/2); 
    Phise = 2*Phise(1:length(freq))/length(Phise); % Use single-sided spectrum.
    % Limit frequency span to 4th harmonic.
    idx = freq<=4*f0; freq = freq(idx); Phise = Phise(idx);

    Iapp = fft(iapp(ss));
    Iapp = 2*Iapp(1:length(freq))/length(Iapp);
    Iapp = Iapp(idx);

    % Extract 1st, 2nd, and 3rd harmonic components.
    fh = (1:3)*f0;
    [~,ih] = min(abs(fh-freq'));
    Phiseh = Phise(ih);
    Iapph = Iapp(ih);
    II = Iapph(1).^(1:3).';
    Phiseh = Phiseh./II;  % Normalize to current amplitude.

    % Store output.
    data.time = time(ss);
    data.time = data.time - data.time(1);
    data.iapp = iapp(ss);
    data.phise = phise(ss);
    data.freq = freq;
    data.Phise = Phise;
    data.Phiseh = Phiseh;
    data.Iapph = Iapph;
end


function RHS = getRHS(t,y,f0,i0,beta,Cdl)
    iapp = I*sin(2*pi*f0*t);
    iF = i0*(exp((1-beta)*f*y)-exp(-beta*f*y));
    RHS = iapp/Cdl - iF/Cdl;
end

end