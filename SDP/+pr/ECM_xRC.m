classdef ECM_xRC < handle
    %ECM_xRC Equivalent-circuit modeling a pulse response waveform.
    % Supports arbitrary series RC circuits.

    properties
        % Circuit Parameters ----------------------------------------------
        % Each of the folling appear in series.

        % Pulse resistance.
        R0 

        % Parallel RC circuits. (taux = 1/RxCx). Vector parameters.
        % For each x, specify any two of (taux,Rx,Cx) for a complete model.
        % Omitted elements are set to NaN.
        nRC     % Number of parallel RC circuits.
        taux
        Rx 
        Cx

        % Total charge-transfer resistance (sum of the Rx).
        Rct

        % Structure of 1sigma error bounds on parameters (if available).
        bound
    end

    methods
        function obj = ECM_xRC(R0, taux, Rx, Cx, varargin)
            p = inputParser;
            p.addOptional('bound',struct,@isstruct);
            p.parse(varargin{:});
            obj.bound = p.Results.bound;

            taux = taux(:);  nTau = length(taux);
            Rx = Rx(:);      nR = length(Rx);
            Cx = Cx(:);      nC = length(Cx);
            nRC = max([nTau,nR,nC]);
            if nTau < nRC, taux = [taux; nan(nRC-nTau,1)]; end
            if nR < nRC, Rx = [Rx; nan(nRC-nR,1)]; end
            if nC < nRC, Cx = [Cx; nan(nRC-nC,1)]; end

            emptytau = isempty(taux) | isnan(taux);
            emptyR = isempty(Rx) | isnan(Rx);
            emptyC = isempty(Cx) | isnan(Cx);

            taux(emptytau) = Rx(emptytau).*Cx(emptytau);
            Rx(emptyR) = taux(emptyR)./Cx(emptyR);
            Cx(emptyC) = taux(emptyC)./Rx(emptyC);

            obj.R0 = R0;
            obj.nRC = nRC;
            obj.taux = taux;
            obj.Rx = Rx;
            obj.Cx = Cx;
            obj.Rct = sum(Rx);
        end

        function v = getV(obj, current, time, components)
            %GETV Compute model-predicted voltage given current and time
            % vectors.

            time = time(:);
            current = current(:);

            if ~exist('component','var')
                all = true;
            else
                components = split(components);
                all = false;
            end

            v = zeros(size(time));
            if all || ismember('R0',components)
                v = v + obj.R0*current;
            end
            if all || ismember('RC',components)
                iRx = pr.ECM_xRC.getiRx(obj.taux,current,time);
                v = v + iRx*obj.Rx;
            end
        end

        function Z = getZ(obj, freq)
            % Compute the impedance of the ECM over the given frequencies.
            Z = ones(size(freq))*obj.R0;
            for k = 1:obj.nRC
                Rk = obj.Rx(k);
                tauk = obj.taux(k);
                Z = Z + Rk./(1 + 1j*2*pi*freq*tauk);
            end
        end
    end

    methods(Static)
        function iRx = getiRx(taux, current, time)
            %GETIRX Compute model-predicted current in Rx (part of the
            % parallel RC pairs). This requires knowlegde of the time
            % constants taux=RxCx only. Multiply each column of iRx by 
            % the corresponding Rx to get the voltage drop across the 
            % corresponding RC pair.
            %
            % NOTE: The time vector must be uniformly sampled! If not, use 
            % interp1 to get uniform spacing.

            % Form discrete-time system using zero-order hold method.
            dt = time(2)-time(1); Ad = exp(-dt./taux(:)'); Bd = 1-Ad;

            % Evalulate output numerically (fastest).
            ik = zeros(1,length(taux)); % Assume iRx(0)=0
            iRx = zeros(length(time),length(taux));
            for k = 1:length(time)
                ik = Ad.*ik + Bd.*current(k);
                iRx(k,:) = ik;
            end
        end

        function vect = pack(model)
            %PACK Stuff the equivalent circuit model instance MODEL into a
            % vector for use in optimization routines. Saves only the time
            % constants.

            vect = model.taux(:);
        end

        function model = unpack(vect, addl)
            %UNPACK Unstuff an equivalent circuit model instance given a
            % vector of time constants VECT, and a structure giving the values of
            % any additional fields, ADDL. Fields that do not appear in either VECT
            % or ADDL are set to nan. 

            if isfield(addl,'R0'), R0 = addl.R0; else, R0 = nan; end
            if isfield(addl,'Rx'), Rx = addl.Rx; else, Rx = nan; end
            if isfield(addl,'Cx'), Cx = addl.Cx; else, Cx = nan; end
            taux = vect(:);

            model = pr.ECM_xRC(R0,taux,Rx,Cx);
        end

        function est = guess(pulse, tmax, nRC)
            %GUESS Generate a rough estimate of model parameters suitable
            %  for initializing a nonlinear optimizer.
            %
            % est = GUESS(pulse, nRC) roughly estimates equivalent-circuit model
            %   parameters from pulse-response data collected in the
            %   laboratory. PULSE is an instance of Pulse containing the
            %   laboratory data. EST contains the parameter estimates (an instance of
            %   EquivCircuit). NRC is the number of parallel RC circuits to fit.

            time = pulse.time;
            voltage = pulse.voltage;
            ipulse = -pulse.ipulse;

            int1 = time<pulse.t0;   % Pre-pulse OCV (OCV1) interval.
            ocv1 = mean(voltage(int1));       % OCV before pulse.

            idx = pulse.t0<=time&time<=tmax;
            time = time(idx); voltage = voltage(idx);
            intss = 0.8*tmax<=time&time<=tmax;       % Steady-state interval (immediately before pulse release).
            intp = pulse.t0<time&time<pulse.tf;      % Pulse interval (instants when current is nonzero).
        
            % Subtract out pre-pulse voltage.
            deltavcal = voltage - ocv1;
            deltavss = mean(deltavcal(intss));
            
            % Estimate R0 and avg(Rx).
            [~, idx0] = min(abs(time));
            R0hat = deltavcal(idx0+1)/ipulse;
            Rbar = (deltavss/ipulse - R0hat)/nRC;
            Rxhat = ones(nRC,1)*Rbar;
            
            % Estimate time constants and Cx.
            deltavcal2 = movmean(deltavcal(intp)-R0hat*ipulse,5);
            deltavcal2 = deltavcal2(:);
            v1tau = ipulse*Rbar*((1:nRC) - exp(-1));
            [~, idx1tau] = min(abs(deltavcal2-v1tau));
            tauhat = time(idx1tau+idx0);
            Cxhat = tauhat./Rxhat;     

            est = pr.ECM_xRC(R0hat, tauhat, Rxhat, Cxhat);
        end

        function [est, J, tp, ip, vp, costFcn] = fit(pulse, guess, tmax, eps, taumin, varargin)
            %FIT Regress an equivalent circuit model to pulse response data 
            % collected in the laboratory.
            %
            % Input:
            %   weighting = array of structures defining weights for the
            %     cost function over time intervals. Should have the
            %     following fields:
            %        t1, t2 - define the time interval over which weight
            %                 applies.
            %        w      - amount by which to weight residuals over the
            %                 time interval (t1, t2). Actual weight is w^2.
            %     The default weight is 1 for unspecified time intervals.

            fieldsOpt = {'taux'};

            % Defaults for positional arguments.
            if ~exist('eps','var')
                eps = 0.5;
            end
            if ~isstruct(eps)
                e = cell(length(fieldsOpt),1);
                e(:) = eps;
                eps = cell2struct(e,fieldsOpt,2);
            end

            % Optional keyword arguments.
            optOptsDefault = optimoptions(@fmincon,'Display','off',...
                    'MaxFunEvals',1e6,'MaxIter',1e3,...
                    'TolFun',1e-20,'TolX',1e-20,'TolCon',1e-20);
            p = inputParser;
            p.addOptional('lb',struct,@isstruct);
            p.addOptional('ub',struct,@isstruct);
            p.addOptional('weighting',[]);
            p.addOptional('optimizationOptions',optOptsDefault);
            p.parse(varargin{:});
            lbUser = p.Results.lb;
            ubUser = p.Results.ub;
            weighting = p.Results.weighting;
            optOpts = p.Results.optimizationOptions;

            % Define optimization bounds. Ensure all parameters >0.
            for k = 1:length(fieldsOpt)
                field = fieldsOpt{k};
                if isfield(lbUser,field)
                    lb.(field) = max(0,lbUser.(field));
                else
                    lb.(field) = max(0,guess.(field).*(1-eps.(field)));
                end
                if isfield(ubUser,field)
                    ub.(field) = max(0,ubUser.(field));
                else
                    ub.(field) = max(0,guess.(field).*(1+eps.(field)));
                end
            end

            % Convert bounds to vectors.
            lb = pr.ECM_xRC.pack(lb); 
            ub = pr.ECM_xRC.pack(ub);

            % Define cost function.
            [costFcn, tp, ip, vp] = ...
                pr.ECM_xRC.costFunctionFactory(pulse,tmax,weighting);

            % Define inequality constraints AX<=B.
            nRC = guess.nRC;
            if nRC > 1
                % Constrain tau1 < tau2 < tau3 < ... to limit search space.
                A = toeplitz([1 zeros(1,nRC-2)],[1 -1 zeros(1,nRC-2)]);
                B = zeros(nRC-1,1);

                % Ensure tau1 > taumin.
                A = [-1 zeros(1,nRC-1); A]; B = [-taumin; B];
            else
                % Ensure tau1 > taumin.
                A = -1; B = -taumin;
            end

            % Optimize!
            vect0 = pr.ECM_xRC.pack(guess);
            vect = fmincon(costFcn,vect0,A,B,[],[],lb,ub,[],optOpts);
            [J, est] = costFcn(vect);
        end

        function [costFcn, tp, ip, vp] = costFunctionFactory(pulse, tmax, weighting)
            %COSTFUNCTIONFACTORY Build a cost function for fitting this
            %  circuit to the specified pulse-response data.
            %
            % costFcn = COSTFUNCTIONFACTORY(pulse,tmax,weighting)
            %  creates a cost function for the data contained in the pulse
            %  resonse PULSE with weighting intervals
            %  defined by WEIGHTING (see fit() for the format of
            %  WEIGHTING). TMAX limits the duration of the pulse to which
            %  the model is fit. COSTFCN is a handle to the cost function.
            %
            % [~, tp, ip, vp] = COSTFUNCTIONFACTORY(...) also returns the
            %  pulse response current IP, voltage VP, and time IP vectors
            %  against which model waveforms are costed.
            
            % NOTE: When modifying this function, take care to note we use
            % the OPPOSITE SIGN CONVENTION for current (charge current is
            % positive). No need to worry about changing the sign
            % convention when calling externally.

            time = pulse.time;
            voltage = pulse.voltage;
            current = pulse.current;

            int1 = time<pulse.t0;     % Pre-pulse OCV (OCV1) interval.
            ocv1 = mean(voltage(int1));

            tp = linspace(0,tmax,round(pulse.fs*tmax));
            ip = interp1(time,-current,tp);     % Use charge (not discharge) current as input! (negative sign)
            vp = interp1(time,voltage-ocv1,tp);

            % Define weighting matrix (faster to store diagonal matrix as vector).
            W = zeros(1,length(tp)); check = zeros(1,length(tp));
            for k = 1:length(weighting)
                indicies = and(weighting(k).t1<tp,tp<weighting(k).t2);
                if any(indicies & check)
                    error('Weighting intervals cannot overlap!');
                end
                W = W + weighting(k).w*indicies;
                check = check | indicies;
            end
            W(W==0) = 1;
            W = W(:);

            costFcn = @cost;

            function [J, model] = cost(tau)
                % Remove duplicate time constants, otherwise the
                % least-squares matrix (H matrix) will be rank deficient.
                % Consider two duplicate time constants tau1 = tau2. The
                % corresponding iR1 and iR2 vectors will be equal, meaning
                % we can write:
                %   deltaV = R1*iR1 + R2*iR2 + ... = (R1+R2)iR + ...
                % That is, we can only identify the sum of R1 and R2 using
                % least squares regression. We'll later use R1=R2=(R1+R2)/2
                % for lack of a better choice.
                [tauUnique,~,idxTauUnique] = unique(tau);
                tauUniqueCounts = sum(tau == tauUnique')';

                % Compute current in R of parallel RC circuits.
                iRx = pr.ECM_xRC.getiRx(tauUnique,ip,tp);

                % Regress R0 and Rx values using weighted ordinary least squares.
                % (Technically should be total least squares due to noisy
                % current measurements, but we'll neglect the noise for now.
                % The current measurements are generally much cleaner than 
                % the voltage, which suffers from high quantization error.)
                params = lsqnonneg(W.*[ip(:) iRx],W.*vp(:));  % The {Rj} are nonnegative, use lsqnonneg!
                R0 = params(1); Rx = params(2:end);

                % Compute weighted residuals.
                resid = W.*(vp(:) - (R0*ip(:) + iRx*Rx));

                % Reintroduce duplicates.
                Rx = Rx./tauUniqueCounts;  % We identified sums of identical resistances, divide to get single resistance.
                Rx = Rx(idxTauUnique);

                % Form complete model estimate.
                model = pr.ECM_xRC.unpack(...
                    tau,struct('R0',R0,'Rx',Rx));

                % Compute cost.
                J = sum(resid.^2);
            end
        end
    end
end