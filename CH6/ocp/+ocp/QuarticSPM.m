classdef QuarticSPM < handle
    %QUARTICSPM Single-particle model with quartic stiochiometry versus r,
    %  lithium-metal chemistry.
    %
    % Note that LMB cells do not have an intercalation negative electrode,
    % and thus only one particle appears in the SPM cell model.

    properties(Constant)
        R = 8.3144598;      % Molar gas constant [J/mol K]
        F = 96485.3329;     % Faraday constant [C/mol]

        % Parameters that are allowed to be optimized.
        FIELDSOPT = {'b', 'm', 'Dref', 'R0', 'R1'};

        % Normalized particle radius at which to evalulate solid
        % diffusivity.
        R_REF = (0.98)^(1/3);
    end

    properties(SetAccess=protected)
        Uocp     % OCP vs stiochiometry function [V].
        dUocp    % Differential capacity vs stiochiometry function.
        zmin     % Lithiation range over which Uocp and dUocp are valid.
        zmax
        T        % Absolute temperature at which OCP data was collected [K].
        b        % Positive-electrode OCP stretch factor.
        m        % Positive-electrode OCP stretch/compress factor.
        Qtilde   % Normalized capacity, Q(tilde) [Ah].
        Dref     % Reference diffusivity (lumped model).
        R0       % Linear resistance parameter [Ohm].
        R1
        Rprofile
        dt       % Sampling interval [sec].

        % Constant state-space matricies.
        Bd
        Cd
    end

    methods
        function obj = QuarticSPM(varargin)
            %QUARTICSPM Construct an instance of the quartic SPM.

            parser = inputParser;
            addParameter(parser,'Uocp',NaN,@(x)isa(x,'function_handle'));
            addParameter(parser,'dUocp',NaN,@(x)isa(x,'function_handle'));
            addParameter(parser,'zmin',NaN,@(x)0<=x&&x<=1);
            addParameter(parser,'zmax',NaN,@(x)0<=x&&x<=1);
            addParameter(parser,'T',NaN,@(x)x>=0);
            addParameter(parser,'b',NaN,@(x)x>=0);
            addParameter(parser,'m',NaN,@(x)x>=0);
            addParameter(parser,'Qtilde',NaN,@(x)x>=0);
            addParameter(parser,'Dref',NaN,@(x)x>=0);
            addParameter(parser,'R0',NaN,@(x)x>=0);
            addParameter(parser,'R1',NaN,@(x)x>=0);
            addParameter(parser,'Rprofile',NaN,@(x)isa(x,'function_handle'));
            addParameter(parser,'dt',NaN,@(x)x>=0);
            parse(parser,varargin{:});

            params = struct2cell(parser.Results);
            idxnan = cellfun(...
                @(x)(~isa(x,'function_handle') && any(isnan(x))),params);
            if any(idxnan)
                error('Missing required parameters %s',...
                      strjoin(params(idxnan),', '));
            end
            
            obj.Uocp = parser.Results.Uocp;
            obj.dUocp = parser.Results.dUocp;
            obj.zmin = parser.Results.zmin;
            obj.zmax = parser.Results.zmax;
            obj.T = parser.Results.T;
            obj.b = parser.Results.b;
            obj.m = parser.Results.m;
            obj.Qtilde = parser.Results.Qtilde;
            obj.Dref = parser.Results.Dref;
            obj.R0 = parser.Results.R0;
            obj.R1 = parser.Results.R1;
            obj.Rprofile = parser.Results.Rprofile;
            obj.dt = parser.Results.dt;

            % Constant discrete-time state-space matricies
            % (Ad, Bd, and Dd are dynamic and computed during simulation).
            obj.Cd = [1, 8/35]; 
        end

        function [voltage, zssk, xk, Dsk] = sim(obj, iapp, soc0, useConstDs)
            %SIM Simulate the SPM, returning the voltage estimate.
            %
            % [voltage,zssk] = SIM(current,soc0) simulates the SPM starting
            %   from the state-of-charge SOC0 given the current vector 
            %   CURRENT, returning the voltage estimate at each sample
            %   in the vector VOLTAGE.
            %
            % The current vector should be sampled at the interval
            % specified in this QuarticSPM object for correct results.

            if ~exist('useConstDs','var')
                useConstDs = false;
            end

            ifk = -iapp;                 % Faradiac current at positive electrode particle surface.
            zssk = zeros(size(iapp));    % Relative stiochiometry at particle surface.
            xk = zeros(2,length(iapp));
            Dsk = zeros(size(iapp));
            voltage = zeros(size(iapp)); % Voltage estimate.

            % Initialize state-space simulation.
            % (State vector consists of average lithiation and average 
            % lithiation gradient.)
            xk(:,1) = [1-soc0; 0];        
            zssk(1) = 1-soc0; 
            voltage(1) = obj.Uocp(zssk(1)) - obj.R0*iapp(1);

            if useConstDs
                z0 = linspace(0,1,1000);
                dUavg = mean((1-obj.m*(obj.b-z0)/obj.b).*(z0-obj.b).*obj.dUocp(z0));
                Ds = abs(dUavg*obj.Dref*obj.F/obj.R/obj.T);
            end

            for k = 1:length(iapp)-1
                % Compute diffusivity for this iteration.
                if useConstDs
                    Dsk(k) = Ds;
                else
                    Rref = obj.R_REF;
                    zss = zssk(k); zavg = xk(1,k); dzavg = xk(2,k);
                    zref = -35*zavg*(1/4 -Rref^2 +(3/4)*Rref^4)...
                        + zss*(39/4 - 35*Rref^2 + (105/4)*Rref^4)...
                        - dzavg*(3 - 10*Rref^2 + 7*Rref^4);
                    Dsk(k) = obj.Dref*obj.F*(1-obj.m*(obj.b-zref)/obj.b)...
                             *(zref-obj.b)*obj.dUocp(zref)/obj.R/obj.T;
                    % Sign varies depending on the sign convention for
                    % dUocp/dz; ensure positive diffusivity.
                    Dsk(k) = abs(Dsk(k));
                end

                % Compute time-varying A(k) and D(k+1) matricies.
                ad = exp(-30*Dsk(k)*obj.dt); 
                bd = (ad-1)/30/Dsk(k)/480/obj.Qtilde;
                Adk = [1 0; 0 ad]; 
                Bdk = [-obj.dt/3600/obj.Qtilde; bd];
                Ddk = -1/Dsk(k)/378000/obj.Qtilde;

                % Propagate state (k => k+1).
                xk(:,k+1) = Adk*xk(:,k) + Bdk*ifk(k);

                % Compute state-space model output (solid surface stiochiometry).
                zssk(k+1) = obj.Cd*xk(:,k+1) + Ddk*ifk(k+1);

                % Saturate lithiation rather than going out of the [0,1] 
                % relative lithiation range. 
                % (We do not have reliable OCP data outside of this range).
                if zssk(k+1) > obj.zmax
                    zssk(k+1) = obj.zmax;
                elseif zssk(k+1) < obj.zmin
                    zssk(k+1) = obj.zmin;
                end

                % Compute voltage estimate.
                voltage(k+1) = obj.Uocp(zssk(k+1)) - obj.R0*iapp(k+1) ...
                    - obj.R1*obj.Rprofile(zssk(k+1))*iapp(k+1);
            end
        end % sim()
    end

    methods(Static)
        function vect = pack(obj, params)
            %PACK Place the specified object fields in a vector.
            %
            %  vect = PACK(obj,params) places specified fields of OBJ into
            %    column vector VECT. PARAMS is a structure whose fieldnames
            %    specify the object fields to stuff into VECT and whose 
            %    values specify (optional) "adapters" (functions) to apply
            %    to each value before placing it into the vector. Each
            %    adapter should be a structure containing the 'pack' field
            %    which specifies the function to apply before packing the
            %    corresponding object field; otherwise, the adapter is
            %    ignored.

            names = fieldnames(params);
            cnt = length(names);

            extra = setdiff(names,QuarticSPM.FIELDSOPT);
            if ~isempty(extra)
                error('Cannot pack parameters %s',strjoin(extra));
            end

            vect = zeros(cnt,1);
            for k = 1:cnt
                name = names{k}; adapter = params.(name);
                if isa(adapter,'struct') && isfield(adapter,'pack')
                    vect(k) = adapter.pack(obj.(name));
                else
                    vect(k) = obj.(name);
                end
            end
        end

        function obj = unpack(obj, vect, params)
            %UNPACK Restore the specified fields from a vector.

            names = fieldnames(params);
            cnt = length(names);

            extra = setdiff(names,QuarticSPM.FIELDSOPT);
            if ~isempty(extra)
                error('Cannot unpack parameters %s',strjoin(extra));
            end

            for k = 1:cnt
                name = names{k}; adapter = params.(name);
                if isa(adapter,'struct') && isfield(adapter,'unpack')
                    obj.(name) = adapter.unpack(vect(k));
                else
                    obj.(name) = vect(k);
                end
            end
        end

        function [spm, data] = fit(spmtest, ocp, QAh, varargin)
            %FIT Regress a SPM to constant-current dis/charge data.
            %
            % [spm, data] = FIT(spmtest,ocp,...) regresses a
            %  single-particle model (SPM) to the discharge data contained
            %  in SPMTEST (an instance of OcpTest) given an estimate of the
            %  cell OCP (either relative or absolute).

            vmin = ocp.vmin;
            vmax = ocp.vmax;
            zmaxU = max(ocp.ZU); zmaxdU = max(ocp.ZdU);
            zminU = min(ocp.ZU); zmindU = min(ocp.ZdU);
            b = zmaxU; % Assume theta=1 when theta(tilde)=1.

            % Specify the parameters we'll optimize below. The diffusivity
            % Dref is typically very small, so we optimize over log10(Dref);
            % the following defines an "adapter" to perform the conversion.
            adapterDref = struct('pack',@(x)log10(x),'unpack',@(x)10.^x);
            paramsopt = struct('Dref',adapterDref,'m',1);
            paramcnt = length(fieldnames(paramsopt));

            parser = inputParser;
            addOptional(parser,'optimizeZMIN',struct.empty,@(x)all(isfield(x,{'lb','ub'})));
            addOptional(parser,'optimizeDREF',struct.empty,@(x)all(isfield(x,{'lb','ub'})));
            addOptional(parser,'vcutoff',vmax-0.5*(vmax-vmin),@(x)x<vmax);
            addOptional(parser,'dt',1,@(x)x>0);
            addOptional(parser,'swarmIterations',25,@(x)x>0);
            addOptional(parser,'swarmSize',50,@(x)x>0);
            addOptional(parser,'fminconIterations',100,@(x)x>0);
            addOptional(parser,'verbose',true,@islogical);
            addOptional(parser,'Rprofile',NaN,@(x)isa(x,'function_handle'));
            parse(parser,varargin{:});

            optZMIN = parser.Results.optimizeZMIN;
            optDREF = parser.Results.optimizeDREF;
            vcutoff = parser.Results.vcutoff;
            dt = parser.Results.dt;
            iterSwarm = parser.Results.swarmIterations;
            pcntSwarm = parser.Results.swarmSize;
            iterFmincon = parser.Results.fminconIterations;
            verbose = parser.Results.verbose;
            Rprofile = parser.Results.Rprofile;
            if isempty(optZMIN)
                error(['Missing optimization bounds for theta(min). ' ...
                       'Provide the argument ''optimizeZMIN''.']);
            end
            if isempty(optDREF)
                error(['Missing optimization bounds for Dref. ' ...
                       'Provide the argument ''optimizeDREF''.']);
            end
            if ~isa(Rprofile,'function_handle')
                error('Missing resistance versus relative lithiation profile');
            end

            % Define optimization bounds.
            zt0 = ocp.zmin;
            z1 = optZMIN.lb; z2 = optZMIN.ub;
            lb.m = (1-z2)/(1-zt0/b); ub.m = (1-z1)/(1-zt0/b);
            lb.Dref = optDREF.lb;    ub.Dref = optDREF.ub;

            % Prepare initial SPM.
            init.T = ocp.temp+273.15;
            init.Qtilde = QAh/(zmaxU-zminU);
            init.Uocp = @ocp.Uocp;
            init.dUocp = @ocp.dUocp;
            init.zmin = min(zminU,zmindU);
            init.zmax = max(zmaxU,zmaxdU);
            init.b = b;
            init.m = (lb.m+ub.m)/2;
            init.Dref = sqrt(lb.Dref*ub.Dref);
            init.R0 = 0; % We compute R0 by linear regression inside the cost function.
            init.R1 = 0;
            init.Rprofile = Rprofile;
            init.dt = dt;
            spmGbl = QuarticSPM(init);

            % Convert optimization bounds to vectors.
            lb = QuarticSPM.pack(lb,paramsopt);
            ub = QuarticSPM.pack(ub,paramsopt);

            % Prepare CC discharge data vectors.
            time = spmtest.script1.time;
            iapp = spmtest.script1.current;
            vtrue = spmtest.script1.voltage;
            disAh = spmtest.script1.disAh;
            dis = iapp > 0;
            time = time(dis); iapp = iapp(dis); vtrue = vtrue(dis); disAh = disAh(dis);

            % SPM does not model terminal voltage well at low SOC;
            % do not regress low-SOC data.
            idx =  vtrue>=vcutoff; %disAh<=0.6*QAh;
            time = time(idx); iapp = iapp(idx); vtrue = vtrue(idx);
            time = time - time(1);

            % Make time vector uniform.
            timeunif = 0:dt:max(time);
            iappGbl = interp1(time,iapp,timeunif,'linear','extrap');
            vtrueGbl = interp1(time,vtrue,timeunif,'linear','extrap');
            timeGbl = timeunif;

            % Determine the starting SOC for the discharge test more 
            % precisely (ideally should be 100%).
            vinitial = spmtest.script1.voltage(spmtest.script1.step==1);
            soc0 = 1-ocp.invUocp(mean(vinitial));

            % Prepare particle swarm matrix. Each row specifies the 
            % parameters associated with a particle in the swarm.
            swarm = zeros(pcntSwarm,paramcnt);
            swarm(1,:) = QuarticSPM.pack(init,paramsopt);
            for k = 1:paramcnt
                swarm(2:end,k) = (ub(k)-lb(k))*rand(pcntSwarm-1,1)+lb(k);
            end

            % Configure particle swarm optimization.
            if verbose
                display = 'iter';
            else
                display = 'off';
            end
            optionsHybrid = optimoptions(@fmincon,'Display',display,...
                'MaxFunEvals',1e6,'MaxIter',iterFmincon,...
                'TolFun',1e-10,'TolX',1e-10,'TolCon',1e-10);
            options = optimoptions(@particleswarm,...
                'Display',display,'UseParallel',false,...
                'FunctionTolerance',1e-10,'SwarmSize',pcntSwarm,...
                'MaxIterations',iterSwarm,'MaxStallIterations',100,...
                'InitialSwarmMatrix',swarm,'FunValCheck','off',...
                'HybridFcn',{@fmincon,optionsHybrid});

            % Run optimization.
            minJ = Inf;
            vect = particleswarm(@cost,paramcnt,lb,ub,options);

            % Output relevant data.
            [J, init.R0, init.R1] = cost(vect);
            init = QuarticSPM.unpack(init,vect,paramsopt);
            spm = QuarticSPM(init);
            data.time = timeGbl;
            data.iapp = iappGbl;
            data.vtrue = vtrueGbl;
            data.soc0 = soc0;
            data.ocp = ocp;
            data.spmtest = spmtest;
            data.J = J;
            
            function [J, R0, R1] = cost(vect)
                 QuarticSPM.unpack(spmGbl,vect,paramsopt);
                 [vest, zss] = spmGbl.sim(iappGbl,soc0);

                 if any(isnan(vest)) || any(isinf(vest))
                    J = Inf;
                    return
                 end

                 vlin = vtrueGbl - vest;
                 Rp = Rprofile(zss(:));
                 Rx = lsqnonneg(-[iappGbl(:),iappGbl(:).*Rp],vlin(:)); % Need to ensure R0>=0!
                 R0 = Rx(1); R1 = Rx(2);
                 vest = vest - iappGbl*R0 - iappGbl.*Rp'*R1;
                 J = rms(vtrueGbl - vest);

                 if verbose && J < minJ
                     minJ = J;
                     fprintf(['J=%10.4f mV, R0=%10.4f mOhm, Dref=%.5e, ' ...
                              'm=%10.4f, b=%10.4f\n'], ...
                              J*1000,R0*1000,spmGbl.Dref,spmGbl.m,spmGbl.b);
                     figure(1); clf;
                     plot(timeGbl,vtrueGbl); hold on;
                     plot(timeGbl,vest);
                     thesisFormat([0.1 0.1 0.5 0.1]);
                     figure(2); clf;
                     plot(timeGbl,vtrueGbl-vest); hold on;
                     yline(0,'k--');
                     thesisFormat([0.1 0.1 0.5 0.1]);
                     drawnow;
                 end
            end
        end % fit()
    end
end