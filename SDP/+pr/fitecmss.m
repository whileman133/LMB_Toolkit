function [params, bounds, sysfit, data] = fitecmss(pulses,varargin)
    %FITECMSS Fit an Equivalent Circuit Model (ECM) to a collection of
    % step-response measurements to estimate the high-frequency 
    % pulse-resistance R0. Uses state-space estimation capability of the
    % system ID toolbox.
    %
    % PARAMS = fitecmss(PULSES) 

    p = inputParser;
    p.addOptional('nRC',2,@isscalar);
    p.addOptional('tmax',100e-6,@isscalar);
    p.addOptional('taumin',1e-6,@isscalar);
    p.parse(varargin{:});
    nRC = p.Results.nRC;
    tmax = p.Results.tmax;
    taumin = p.Results.taumin;

    % Collect pulse data.
    Udata = cell(1,length(pulses));
    Ydata = cell(1,length(pulses));
    RLargeSignal = zeros(1,length(pulses));
    kk = zeros(1,length(pulses));
    tp = linspace(0,tmax,round(pulses(1).fs*tmax));
    Ts = tp(2)-tp(1);
    for k = 1:length(pulses)
        pulse = pulses(k);
        
        time = pulse.time;
        voltage = pulse.voltage;
        current = pulse.current;

        int1 = time<pulse.t0;   % Pre-pulse OCV (OCV1) interval.
        ocv1 = mean(voltage(int1));

        % Use charge (not discharge) current as input! (negative sign)
        ip = interp1(time,-current,tp); 
        vp = interp1(time,voltage-ocv1,tp);

        Udata{k} = ip(:);
        Ydata{k} = vp(:)/Ts;
        RLargeSignal(k) = vp(end)/ip(end);
        sigmaV = std(voltage(int1));
        sigmaI = std(current(int1));
        kk(k) = Ts*sigmaI/sigmaV;
    end
    pulsedata = iddata(Ydata,Udata,Ts);

    % Create initial guess.
    ecminitial = pr.ECM_xRC.guess(pulses(1),tmax,nRC);
    taux0 = ecminitial.taux;
    Rx0 = ecminitial.Rx;
    R0 = ecminitial.R0;

    % Construct initial state-space model.
    A = diag(-1./taux0);
    B = ones(nRC,1);
    C = (Rx0(:)./taux0(:))'/Ts;
    D = R0/Ts;
    K = ones(nRC,1)*mean(kk);
    x0 = zeros(nRC,1);
    sysinitial = idss(A,B,C,D,K,x0,0);

    % Define optimization bounds.
    sysinitial.Structure.A.Free = (A~=0);
    sysinitial.Structure.B.Free = false;
    sysinitial.Structure.C.Free = true;
    sysinitial.Structure.D.Free = true;
    sysinitial.Structure.K.Free = false;
    for k = 1:nRC
        sysinitial.Structure.A.Minimum(k,k) = -1/taumin;
        sysinitial.Structure.A.Maximum(k,k) = 0;
        sysinitial.Structure.C.Minimum(1,k) = 0/Ts;
        sysinitial.Structure.C.Maximum(1,k) = max(RLargeSignal)/taumin/Ts;
    end
    sysinitial.Structure.D.Minimum = 0/Ts;
    sysinitial.Structure.D.Maximum = max(RLargeSignal)/Ts;

    % Perform model regression.
    opt = ssestOptions('InitialState','zero');
    opt.SearchOptions.Tolerance = 1e-9;
    opt.SearchOptions.MaxIterations = 1e3;
    sysfit = ssest(pulsedata,sysinitial,opt);
    [param, bound] = getpvec(sysfit,'free');

    % Collect regressed parameter values.
    params.tau = -1./param(1:nRC);
    params.Rx = Ts*param(nRC+1:nRC+2).*params.tau;
    params.R0 = Ts*param(nRC+3);
    if ~isempty(bound)
        bounds.R0 = Ts*bound(nRC+3);
    else
        bounds = struct;
        warning('Error bounds not available');
    end

    % Collect pulse waveforms.
    data.Vp = cellfun(@(vp)vp*Ts,Ydata,'UniformOutput',false);
    data.Ip = Udata;
    data.tp = tp;
    data.Ts = Ts;
end