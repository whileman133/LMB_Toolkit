classdef MSMR < handle
    %MSMR Multiple-site multiple-reaction (MSMR) OCP model.

    properties(Constant)
        R = 8.3144598;      % Molar gas constant [J/mol K]
        F = 96485.3329;     % Faraday constant [C/mol]
    end

    properties(SetAccess=protected)
        J        % Number of galleries.
        Xj       % Fraction of sites belonging to each gallery.
        Uj0      % Standard electrode potentials for each reaction.
        Wj       % Nernst equation shape factor for each gallery.
        zmin     % Minimum lithiation (could be absolute or relative depending on the context).
        zmax     % Maximum lithiation (could be absolute or relative depending on the context).
        name     % Model name (if applicable).
        ref      % Citation/reference (if applicable).
    end

    methods
        function obj = MSMR(Xj, Uj0, Wj, zmin, zmax, name, reference, sortParams)
            %MSMR Construct an MSMR object with the specified parameters.

            if ~exist('name','var')
                name = '';
            end
            if ~exist('reference','var')
                reference = '';
            end
            if ~exist('sortParams','var')
                sortParams = true;
            end

            obj.Xj = Xj(:);
            obj.Uj0 = Uj0(:);
            obj.Wj = Wj(:);
            obj.J = length(Xj);
            obj.zmin = zmin;
            obj.zmax = zmax;
            obj.name = name;
            obj.ref = reference;

            % Sort parameters by equilibrium potential.
            if sortParams
                [obj.Uj0, idx] = sort(obj.Uj0,'descend');
                obj.Xj = obj.Xj(idx);
                obj.Wj = obj.Wj(idx);
            end
        end

        function ocp = asOCP(obj, vmin, vmax, T)
            %ASOCP Convert this MSMR model to a standard OCP object.

            [U, dzdv, Z] = obj.getOCP('volt', vmin, vmax, 1000, T);
            z1 = min(Z);
            z2 = max(Z);
            Z = (Z - z1)/(z2-z1);
            ocp = OCP(vmin, vmax, T, Z, U, Z, 1./dzdv, 'MSMR');
        end

        function modNew = delete(obj, galleries)
            %DELETE Remove galleries from a model, returing a new model.
            %
            % modNew = DELETE(galleries) removes the galleries with
            %  indicies GALLERIES and returns the resulting model.
            %  The the Xj are updated to ensure sum(Xj)=1.

            j = 1:obj.J;
            idx = any(j(:)==galleries(:)',2);

            XjNew = obj.Xj(~idx);
            Uj0New = obj.Uj0(~idx);
            WjNew = obj.Wj(~idx);

            resid = 1 - sum(XjNew);
            XjNew = XjNew + resid/length(XjNew);
            modNew = ocp.MSMR(XjNew, Uj0New, WjNew, obj.zmin, obj.zmax);
        end

        function modNew = insert(obj, Uj0p, varargin)
            %INSERT Insert galleries into a model, returning a new model.

            p = inputParser;
            p.addOptional('Xj',0.1*ones(size(Uj0p)),@(x)all(0<=x&x<=1));
            p.addOptional('Wj',ones(size(Uj0p)));
            p.parse(varargin{:});
            Xjp = p.Results.Xj;
            Wjp = p.Results.Wj;

            Uj0New = [obj.Uj0; Uj0p(:)];
            XjNew = [obj.Xj-sum(Xjp)/obj.J; Xjp(:)];
            WjNew = [obj.Uj0; Wjp(:)];
            modNew = ocp.MSMR(XjNew, Uj0New, WjNew, obj.zmin, obj.zmax);
        end

        function [U, theta, xj, T, npoints, dUdtheta] = ocp(obj, varargin)
            %OCP Calculate OCP predicted by the MSMR model.
            %
            % [U, theta] = OCP() computes the MSMR OCP curve over the
            %   default lithiation range (zmin to zmax).
            %   Potential points are evenly spaced.
            %
            % [U, theta] = OCP('vmin',vmin,'vmax',vmax) computes the MSMR
            %   OCP curve over the potential range VMIN to VMAX.
            %   Potential points are evenly spaced.
            %
            % [U, theta] = OCP('thetamin',tmin,'thetamax',tmax) computes 
            %   the MSMR OCP curve over the lithiation range TMIN to TMAX.
            %   Potential points are evenly spaced. Lithiation points are
            %   *unevely* spaced.
            %
            % [U, theta] = OCP('voltage',voltage) computes the MSMR
            %   OCP curve at potential points given by the vector VOLTAGE.
            %
            % [U, theta] = OCP('theta',theta) computes the MSMR
            %   OCP curve at lithiation points given by the vector THETA.
            %   NPOINTS determines the accuracy of the approximation.
            %
            % [...] = OCP(...,'T',T) computes the OCP at temperature T degC.
            %
            % [...] = OCP(...,'npoints',npoints) uses NPOINTS points to
            %   generate the potential vector (if not supplied).
            %
            % [..., xj] = OCP(...) also returns the partial lithiations of
            %   the galleries, an intermediate computation. XJ is a matrix
            %   whose rows correspond to galleries and whose columns
            %   correspond to individual points.

            p = inputParser;
            p.addOptional('vmin',[]);
            p.addOptional('vmax',[]);
            p.addOptional('voltage',[]);
            p.addOptional('thetamin',[]);
            p.addOptional('thetamax',[]);
            p.addOptional('theta',[]);
            p.addOptional('T',25);
            p.addOptional('npoints',10000);
            p.parse(varargin{:});
            vmin = p.Results.vmin;
            vmax = p.Results.vmax;
            voltage = p.Results.voltage;
            thetamin = p.Results.thetamin;
            thetamax = p.Results.thetamax;
            theta = p.Results.theta;
            T = p.Results.T;
            npoints = p.Results.npoints;

            % Ensure everything is a row vector.
            voltage = voltage(:)';
            theta = theta(:)';

            % Determine potential vector.
            if ~isempty(vmin) && ~isempty(vmax)
                U = linspace(vmin, vmax, npoints);
            elseif ~isempty(voltage)
                U = voltage;
            elseif ~isempty(thetamin) && ~isempty(thetamax)
                vmax = obj.lith2ocp(thetamin,T);
                vmin = obj.lith2ocp(thetamax,T);
                U = linspace(vmin, vmax, npoints);
            elseif ~isempty(theta)
                vmax = obj.lith2ocp(min(theta),T);
                vmin = obj.lith2ocp(max(theta),T);
                U = linspace(vmin, vmax, npoints);
            else
                vmax = obj.lith2ocp(obj.zmin,T);
                vmin = obj.lith2ocp(obj.zmax,T);
                U = linspace(vmin, vmax, npoints);
            end

            % Compute stoichiometry vector.
            f = obj.f(T);
            expt = exp(f*(U-obj.Uj0)./obj.Wj);
            xj = obj.Xj./(1+expt);
            dxjdU = -f.*obj.Xj.*expt./obj.Wj./(1+expt).^2;
            Z = sum(xj,1);
            dUdtheta = 1./sum(dxjdU,1);

            % Interpolate over requested stiochiometry vector (if supplied).
            if ~isempty(theta)
                Unew = interp1(Z,U,theta,'linear','extrap');
                dUdzNew = interp1(Z,dUdtheta,theta,'linear','extrap');
                xjNew = zeros(obj.J,length(Unew));
                for j = 1:obj.J
                    xjNew(j,:) = interp1(U,xj(j,:),Unew,'linear','extrap');
                end
                U = Unew;
                dUdtheta = dUdzNew;
                xj = xjNew;
            else
                theta = Z;
            end
        end

        function [Rct, i0, i0j, U, theta, xj] = Rct(obj, params, varargin)
            %RCT Compute the MSMR charge-transfer resistance versus SOC.
            %
            % Rct = RCT(params,...) computes the MSMR charge-transfer
            %   resistance at the SOC setpoints specified in (...) using
            %   exactly the same input syntax as MSMR.OCP(). PARAMS is
            %   a structure containing vector fields specifying the 
            %   kinetics parameters k0 and alpha.
            %
            % [..., i0, i0j] = RCT(...) also returns the total exchange
            %   current density i0 and that associated with each gallery i0j.
            %
            % [..., U, theta, xj] = RCT(...) also returns the potential, 
            %   lithiation, and partial lithiation corresponding to the SOC 
            %   setpoints.
            %
            % This function assumes theta_e=1 to perform the computation.
            % (Holds for small-signal linear approximation used in transfer 
            % function model).
            
            % Determine partial lithiations xj.
            [U, theta, xj, T] = obj.ocp(varargin{:});

            % Collect parameters.
            k0 = params.k0;        
            alpha = params.alpha; 
            X = obj.Xj(:);
            omega = obj.Wj(:);

            % If function-based parameters are supplied, evaulate.
            if isa(k0,'function_handle'), k0 = k0(0,T); end
            if isa(alpha,'function_handle'), alpha = alpha(0,T); end
            k0 = k0(:);
            alpha = alpha(:);

            % Now compute the charge-transfer resistance.
            i0j = k0.*xj.^(omega.*alpha).*(X-xj).^(omega.*(1-alpha));
            i0 = sum(i0j,1);
            Rct = 1./i0/obj.f(T);
        end

        function [Uocp, dzdv, Z, xj] = getOCP(obj, method, vzmin, vzmax, npoints, T)
            %GETOCP Assuming equilibrium so that Uj = Uocp, compute partial
            % surface stiochiometries corresponding to each Uocp for
            % each gallery. Also evalulate total surface stoichiometry
            % and differential capacity.
            %
            % DEPRECIATED. Use MSMR.OCP() instead!
            %
            % Input:
            %   vmin, vmax = potential range over which to evalulate OCP
            %   npoints = number of OCP values to evalulate between vmin, vmax
            %   T = temperature, defaults to 25 [degC]
            %
            % Output:
            %   Uocp = open-circuit potential vector.
            %   dzdv = differential capacity vector.
            %   Z = stiochiometry vector.
            %   xj = matrix composed of column stiochiometry vectors for each gallery.

            % Defaults ----------------------------------------------------
            if ~exist('method','var')
                method = 'volt';
            end

            if strcmp(method,'lith')
                if ~exist('vzmin','var')
                    vmin = obj.lith2ocp(obj.zmax,T);
                else
                    vmin = obj.lith2ocp(vzmax,T);
                end
                if ~exist('vzmax','var')
                    vmax = obj.lith2ocp(obj.zmin,T);
                else
                    vmax = obj.lith2ocp(vzmin,T);
                end
            elseif strcmp(method,'volt')
                vmin = vzmin;
                vmax = vzmax;
            else
                error('Choose method=''lith'' or method=''volt''');
            end

            if ~exist('npoints','var')
                npoints = 1000;
            end
            if ~exist('T','var')
                T = 25;
            end

            % Computation -------------------------------------------------
            f = obj.f(T);
            Uocp = linspace(vmin, vmax, npoints);
            expt = exp(f*(Uocp-obj.Uj0)./obj.Wj);
            xj = obj.Xj./(1+expt);
            dxjdU = -f.*obj.Xj.*expt./obj.Wj./(1+expt).^2;
            Z = sum(xj,1);
            dzdv = sum(dxjdU,1);
        end

        function plot(obj, vmin, vmax, npoints)
            %PLOT Plot the OCP and differential capacity over the
            %  specified range.

            [Uocp, dthetadU, theta, xj] = obj.getOCP('volt', vmin, vmax, npoints);

            figure; tiledlayout(3,1);
            nexttile;
            plot(xj, Uocp, '--'); hold on;
            plot(theta, Uocp, 'k');
            xlabel('Surface Stoichiometry, \theta_{se}');
            ylabel('Open-Circuit Potential, U_{ocp} [V]');
            title(sprintf('%s OCP', obj.name));
            nexttile;
            plot(Uocp, abs(dthetadU), 'k');
            xlabel('Open-Circuit Potential, U_{ocp} [V]');
            ylabel('|d\theta_{se}/dV| [V^{-1}]');
            title('Differential Capacity vs OCP');
            nexttile;
            plot(theta, abs(dthetadU), 'k');
            xlabel('Surface Stoichiometry, \theta_{se}');
            ylabel('|d\theta_{se}/dV| [V^{-1}]');
            title('Differential Capacity vs Stoichiometry');
        end

        function [Uocp, xj] = lith2ocp(obj, z, T)
            %LITH2OCP Calculate the OCP corresponding to a
            %  particular surface lithiation by solving the nonlinear
            %  equation in Uocp. This implementation uses fzero.
            %
            % Input:
            %   z = surface lithiation / stiochiometry
            %   T = temperature, defaults to 25 [degC]
            %
            % Output:
            %   Uocp = open-circuit potential corresponding to z
            %   xj =   partial surface stiochiometries corresponding to z
            %          and Uocp

            if ~exist('T','var')
                T = 25;
            end
            f = obj.f(T);

            Uocp = fzero(@lithiationBalanceEquation, mean(obj.Uj0));
            xj = obj.Xj./(1+exp(f*(Uocp-obj.Uj0)./obj.Wj));
            
            function residual = lithiationBalanceEquation(U)
                residual = z - sum(obj.Xj./(1+exp(f*(U-obj.Uj0)./obj.Wj)));
            end
        end

        function print(obj, title)
            %PRINT Print a string representation of the MSMR model to the
            % command window.

            if exist('title','var')
                fprintf('MSMR MODEL: %s\n', title);
            else
                fprintf('MSMR MODEL\n');
            end
            fprintf('\tj =   ');
            fprintf('%-7d', 1:obj.J);
            fprintf('\n');
            fprintf('\tUj0 = '); 
            fprintf('%-7.3f', obj.Uj0(:)');
            fprintf('\n');
            fprintf('\tXj =  ');
            fprintf('%-7.3f', obj.Xj(:)');
            fprintf('\n');
            fprintf('\tWj =  ');
            fprintf('%-7.3f', obj.Wj(:)');
            fprintf('\n');
        end
    end

    methods(Static)
        function vect = pack(model, fixedParams)
            %PACK Reduce an MSMR model to a vector of parameters. For use,
            % for example, with optimization routines.
            %
            % VECT = PACK(MODEL) compacts parameters in MODEL into a 
            %   column vector. MODEL may be an instance of MSMR or a
            %   structure containing fields Xj,Uj0,Wj,zmin,zmax (useful for
            %   making vectors of upper and lower bounds).
            %
            % VECT = PACK(...,FIXEDPARAMS) excludes the fields of the
            %   structure FIXEDPARAMS from the parameter vector.

            if ~exist('fixedParams','var')
                fixedParams = struct;
            end

            if isfield(fixedParams,'Xj'), Xj = []; else, Xj = model.Xj; end
            if isfield(fixedParams,'Uj0'), Uj0 = []; else, Uj0 = model.Uj0; end
            if isfield(fixedParams,'Wj'), Wj = []; else, Wj = model.Wj; end
            if isfield(fixedParams,'zmin'), zmin = []; else, zmin = model.zmin; end
            if isfield(fixedParams,'zmax'), zmax = []; else, zmax = model.zmax; end

            vect = [Xj; Uj0; Wj; zmin; zmax];
        end

        function model = unpack(vect, fixedParams)
            %UNPACK Restore an MSMR model from the given parameter vector.
            % For use, for example, with optimization routines.

            if ~exist('fixedParams','var')
                fixedParams = struct;
            end

            numVects = 3; numScalars = 2;
            if isfield(fixedParams,'Xj'), numVects = numVects - 1; Xj = fixedParams.Xj; J = length(fixedParams.Xj); end
            if isfield(fixedParams,'Uj0'), numVects = numVects - 1; Uj0 = fixedParams.Uj0; J = length(fixedParams.Uj0); end
            if isfield(fixedParams,'Wj'), numVects = numVects - 1; Wj = fixedParams.Wj; J = length(fixedParams.Wj); end
            if isfield(fixedParams,'zmin'), numScalars = numScalars - 1; zmin = fixedParams.zmin; end
            if isfield(fixedParams,'zmax'), numScalars = numScalars - 1; zmax = fixedParams.zmax; end

            cursor = 1;
            if ~exist('J','var'), J = (length(vect)-numScalars)/numVects; end
            if ~exist('Xj','var'), Xj = vect(cursor:J+cursor-1); cursor = cursor + J; end
            if ~exist('Uj0','var'), Uj0 = vect(cursor:J+cursor-1); cursor = cursor + J; end
            if ~exist('Wj','var'), Wj = vect(cursor:J+cursor-1); cursor = cursor + J; end
            if ~exist('zmin','var'), zmin = vect(cursor); cursor = cursor + 1; end
            if ~exist('zmax','var'), zmax = vect(cursor); end

            % sort=false since parameters will already be sorted if
            % unpacking from pack()
            model = ocp.MSMR(Xj,Uj0,Wj,zmin,zmax,'','',false);
        end

        function f = f(T)
            %F Compute f=F/(R*T) for the supplied T in [C].

            if ~exist('T','var')
                T = 25;
            end
            F = ocp.MSMR.F;
            R = ocp.MSMR.R;
            f = F/R/(T+273.15);
        end

        function obj = NMC622()
            %NMC622 Construct a model for NMC622 using Verbrugge, 2017.
            name = "NMC622//Li Electrode";
            reference = ...
                "Mark Verbrugge et al. 2017 J. Electrochem. Soc. 164 E3243";
            obj = ocp.MSMR(...
                [0.13442 0.3246 0.21118 0.32980]', ...
                [3.62274 3.72645 3.90575 4.22955]', ...
                [0.96710 1.39712 3.505 5.52757]', ...
                1, 0.1325, name, reference);
        end

        function ax = compare(umin, umax, dzmax, T, varargin)
            %COMPARE Compare OCP data collected in the laboratory and
            % MSMR models.

            model = ocp.MSMR.empty;       modelNames = cell.empty;
            lab = ocp.Estimate.empty;     labNames = cell.empty;

            if length(varargin) >= 4
                if mod(length(varargin),2) == 1
                    % Odd number of arguments.
                    error('Expected even number of variable arguments.');
                else
                    % Even number of arguments.
                    names = varargin(1:2:length(varargin));
                    objs = varargin(2:2:length(varargin));
                end
            else
                error('Supply at least two objects for comparison.');
            end

            for k = 1:length(objs)
                name = names{k};
                obj = objs{k};

                if isa(obj,'ocp.MSMR')
                    model(end+1) = obj;
                    modelNames{end+1} = name;
                elseif isa(obj,'ocp.Estimate')
                    lab(end+1) = obj;
                    labNames{end+1} = name; 
                else
                    error('Second argument must be an instance of ocp.MSMR or ocp.Estimate!')
                end
            end

            figure;
            ax(1) = subplot(1,2,1);
            ax(2) = subplot(1,2,2);

            if isempty(lab)
                % Compare models.
                for m = model
                    [modV, moddzdv, modZ, ~] = ... 
                        m.getOCP('lith',0.0001,0.9999,1000,T);
                    subplot(1,2,1); plot(modZ, modV); hold on;
                    subplot(1,2,2); plot(modV, abs(moddzdv)); hold on;
                    subplot(1,2,3); plot(modZ, abs(moddzdv)); hold on;
                end
                subplot(1,2,1); legend(modelNames);
            elseif isempty(model)
                % Compare lab measurements.

                % Filter lab data to within voltage limits.
                data(length(lab)) = struct('V',[],'Z',[],'refV',[],'refZ',[],'dzdvRefV',[],'dzdvRefZ',[]);
                for k = 1:length(lab)
                    vmin = lab(k).ocptest.vmin;
                    vmax = lab(k).ocptest.vmax;
                    idxV = and(lab(k).V<=vmax,lab(k).V>=vmin);
                    idxRefV = and(lab(k).refV<=vmax,lab(1).refV>=vmin);
                    data(k).V = lab(k).V(idxV); 
                    data(k).Z = lab(k).Z(idxV);
                    data(k).refV = lab(k).refV(idxRefV);
                    data(k).dzdvRefV = lab(k).dzdvRefV(idxRefV);
                    data(k).refZ = lab(k).refZ;
                    data(k).dzdvRefZ = lab(k).dzdvRefZ;
                end

                % Plot results.
                for k = 1:length(lab)
                    subplot(1,2,1); plot(data(k).Z, data(k).V); hold on;
                    subplot(1,2,2); plot(data(k).refV, abs(data(k).dzdvRefV)); hold on;
                    subplot(1,2,3); plot(data(k).refZ, data(k).dzdvRefZ); hold on;
                end
                subplot(1,2,1); legend(labNames);
            elseif length(model) == 1
                % Compare a model with lab measurements.

                T = lab(1).ocptest.temp;
                data(length(lab)) = struct('V',[],'Z',[],'refV',[],'refZ',[],'dzdvRefV',[],'dzdvRefZ',[]);
                for k = 1:length(lab)
                    % Filter lab data to within voltage limits.
                    vmin = lab(k).ocptest.vmin;
                    vmax = max(lab(k).V); %lab(k).ocptest.vmax;
                    idxV = and(lab(k).V<=vmax,lab(k).V>=vmin);
                    idxRefV = and(lab(k).refV<=vmax,lab(k).refV>=vmin);
                    data(k).V = lab(k).V(idxV); 
                    data(k).Z = lab(k).Z(idxV);
                    data(k).refV = lab(k).refV(idxRefV);
                    data(k).dzdvRefV = lab(k).dzdvRefV(idxRefV);
                    data(k).refZ = lab(k).refZ;
                    data(k).dzdvRefZ = lab(k).dzdvRefZ;

                    % Convert lab data to absolute stiochiometry for comparison
                    % with the model.
                    data(k).Z = model.zmin + data(k).Z*(model.zmax - model.zmin);
                    data(k).refZ = model.zmin + data(k).refZ*(model.zmax - model.zmin);
                    data(k).dzdvRefV = data(k).dzdvRefV*(model.zmax - model.zmin);
                    data(k).dzdvRefZ = data(k).dzdvRefZ*(model.zmax - model.zmin);
                end
    
                [modV, moddzdv, modZ, modxj] = ... 
                    model.getOCP('volt',umin,umax,1000,T);

                % Plot results.
                subplot(1,2,1); plot(modZ, modV); hold on;
                subplot(1,2,2); plot(modV, abs(moddzdv)); hold on;
                for k = 1:length(lab)
                    subplot(1,2,1); plot(data(k).Z, data(k).V);
                    subplot(1,2,2); plot(data(k).refV, abs(data(k).dzdvRefV));
                end
                subplot(1,2,1); plot(modxj, modV, '--');

                % Plot markers indicating the range of the collected data.
                colors = get(gca,'colororder');
                for k = 1:length(lab)
                    subplot(1,2,1); 
                    plot(data(k).Z(1),data(k).V(1),'o','LineWidth',10,'Color',colors(k+1,:)); 
                    plot(data(k).Z(end),data(k).V(end),'o','LineWidth',10,'Color',colors(k+1,:)); 
                    subplot(1,2,2); 
                    plot(data(k).refV(1),data(k).dzdvRefV(1),'o','LineWidth',10,'Color',colors(k+1,:)); 
                    plot(data(k).refV(end),data(k).dzdvRefV(end),'o','LineWidth',10,'Color',colors(k+1,:));
                end

                subplot(1,2,1); legend([modelNames(:)',labNames(:)']);
            else
                error('Unsupported argument mix.')
            end

            subplot(1,2,1);
            xlabel('Stoichiometry, \theta');
            ylabel('Open-Circuit Potential, U_{ocp} [V]');
            title('OCP vs Stiochiometry');
            ylim([umin umax]);
            subplot(1,2,2);
            xlabel('Open-Circuit Potential, U_{ocp} [V]');
            ylabel('|d\theta/dU| [V^{-1}]');
            title('Differential Capacity vs OCP');
            xlim([umin umax]); ylim([0 dzmax]);
        end

        function [Xj, U0j, Wj, dv, err] = guess(ocpest, maxJ, output, dvInitial, dvIncrement, dvFinal)
            %GUESS Attempt to identify MSMR parameters from lab data over
            %  relative lithiation.
            %
            %  The reported Xj values should be scaled by 
            %  [theta_max-theta_min] when converting to absolute lithiation. 
            %  Manually adjusting the Xj is neccesary to ensure sum(Xj)=1 
            %  after scaling. In particular, consider stretching the Xj of 
            %  the galleries at the "edges" of the differential capacity 
            %  curves whose Xj may be larger than shown in collected data.
            %
            % This works well only when all of the galleries are observable
            % in the OCP curve. Manual fudging of this guess is required 
            % for materials such as NMC.
            %
            % Input:
            %   ocpest = OCP estimate object containing the relevant data.
            %   maxJ   = maximum number of galleries allowed (used to adjust
            %            differential capacity smoothing.)
            %   output = 'verbose' to show console outputs and plots
            %            'quiet' to suppress output
            %
            % Output:
            %   Xj, U0j, Wj = vectors of MSMR parameters
            %   dv = the bin size used for the differential capacity
            %     computation
            %   err = true if a problem was encountered guessing any of the
            %     parameters

            % Defaults.
            if ~exist('output','var')
                output = 'quiet';
            end
            if ~exist('dvFinal','var')
                dvFinal = 0.1;
            end
            if ~exist('dvIncrement','var')
                dvIncrement = 1e-3;
            end
            if ~exist('dvInitial','var')
                dvInitial = 1e-3;
            end

            if strcmp(output,'verbose')
                verbose = true;
            else
                verbose = false;
            end

            f = ocp.MSMR.f(ocpest.ocptest.temp);
            idx = and(ocpest.V<=ocpest.ocptest.vmax,ocpest.V>=ocpest.ocptest.vmin);
            theta = ocpest.Z(idx); Uocp = ocpest.V(idx);

            % Compute differential capacity with progressively incresing
            % smoothing until we have less than maxJ galleries.
            for dv = dvInitial:dvIncrement:dvFinal
                [dzdv_V, V, dzdv_Z, Z] = smoothdiff(theta, Uocp, dv);
                [~, idxVPeaks] = findpeaks(dzdv_V);
                [~, idxZValley] = findpeaks(-dzdv_Z);
                if(length(idxVPeaks) <= maxJ && ...
                   length(idxZValley) <= length(idxVPeaks)-1)
                    break;
                end
            end

            % Determine location of peaks over voltage, sort descending.
            vPeaks = V(idxVPeaks);
            dzdvPeaks = dzdv_V(idxVPeaks);
            [vPeaks, vSortIdx] = sort(vPeaks,'descend');
            dzdvPeaks = dzdvPeaks(vSortIdx);

            % Determine locations of "valleys" over stiochometry, sort
            % ascending.
            zValley = sort(Z(idxZValley));

            % Estimate MSMR model parameters.
            U0j = vPeaks;
            Xj = diff([0 zValley 1]);

            if length(Xj) ~= length(U0j)
                Wj = []; 
                err = true;
                warning(['Xj-U0j mismatch. Could not estimate Wj. ' ...
                    'Incomplete MSMR estimate. ' ...
                    'Try changing the number of galleries.']);
            else
                Wj = f*Xj./dzdvPeaks/4;
                err = false;
            end

            if verbose
                if ~err
                    estMSMR = ocp.MSMR(Xj,U0j,Wj,0,1);
                    [estUocp, estdthetadU, estTheta, estxj] = ... 
                        estMSMR.getOCP( ...
                            'volt',ocpest.ocptest.vmin, ocpest.ocptest.vmax, ... 
                             1000, ocpest.ocptest.temp);
                end
    
                fprintf('MSMR PARAMETER GUESS\n');
                fprintf('\tdv = %.3fV\n', dv);
                fprintf('\tU0j = '); 
                fprintf('%.3fV ', U0j);
                fprintf('\n');
                fprintf('\tXj = ');
                fprintf('%.3f ', Xj);
                fprintf('\n');
                fprintf('\tWj = ');
                fprintf('%.3f ', Wj);
                fprintf('\n');
    
                figure;
    
                subplot(1,3,1);
                if ~err
                    plot(estTheta, estUocp); hold on;
                end
                plot(theta, Uocp); hold on;
                if ~err
                    plot(estxj, estUocp, '--');
                end
                xlabel('Stoichiometry, \theta');
                ylabel('Open-Circuit Potential, U_{ocp} [V]');
                title('OCP vs Stiochometry');
                legend('Guess', 'Lab');
                
                subplot(1,3,2);
                if ~err
                    plot(estUocp, abs(estdthetadU)); hold on;
                end
                plot(V, dzdv_V); hold on;
                xline(V(idxVPeaks), '--');
                xlabel('Open-Circuit Potential, U_{ocp} [V]');
                ylabel('|d\theta/dV| [V^{-1}]');
                title('Differential Capacity vs OCP');
    
                subplot(1,3,3);
                if ~err
                    plot(estTheta, abs(estdthetadU)); hold on;
                end
                plot(Z, dzdv_Z); hold on;
                xline(Z(idxZValley), '--');
                xlabel('Stoichiometry, \theta');
                ylabel('|d\theta/dV| [V^{-1}]');
                title('Differential Capacity vs Stoichiometry');
    
                thesisFormat([0.2 0.1 0.1 0.1]);
            end
        end

        function [est, Jcost] = fit(ocpestimates, J, vmin, vmax, varargin)
            %FIT Regress an MSMR model to laboratory OCP measurements.
            %
            % [EST, J] = MSMR.FIT(OCPEST, GUESS, VMIN, VMAX, EPS, W, OPTS) regresses 
            %  an MSMR model to the OCP estimate OCPEST starting from the 
            %  initial guess GUESS (an instance of MSMR) over the voltage 
            %  range defined by VMIN and VMAX. The MSMR parameters are
            %  constrained within the initial guess values by plus or minus
            %  EPS*100 percent. W is the wieght on the agreement of the 
            %  data-derived and model-derived differential capacity curves.
            %  OPTS are the fmincon optimization options (optional).
            %  EST is the best-fit MSMR model found. J is the cost
            %  associated with the best fit.
            %
            % EPS may be a scalar or a structure. A scalar applies uniform
            % boundaries to each parameter. For nonuniform boundaries,
            % provide structure of values with fields Xj,Uj0,Wj,zmin,zmax.
            % Values corresponding to Xj,Uj0,Wj may be either scalars or
            % vectors (for nonuniform boundaries on each component).
            %
            % Use the MSMR.guess() function to generate an initial guess 
            % from lab data.

            msmrFields.Xj.len = J;
            msmrFields.Uj0.len = J;
            msmrFields.Wj.len = J;
            msmrFields.zmin.len = 1;
            msmrFields.zmax.len = 1;
            msmrFieldnames = fieldnames(msmrFields);
            msmrFieldnamesLith = {'Xj','zmin','zmax'};

            % Parse input.
            parser = inputParser;
            parser.addParameter('initial',[],@(x)isa(x,'ocp.MSMR')||isempty(x));
            parser.addParameter('eps',0.5,@(x)isstruct(x)||x>=0);
            parser.addParameter('lb',struct,@isstruct);
            parser.addParameter('ub',struct,@isstruct);
            parser.addParameter('fix',struct,@isstruct);
            parser.addParameter('Usep',0,@(x)x>=0);
            parser.addParameter('w',0.001,@(x)x>=0);
            parser.addParameter('gaPopulationSize',200,@(x)x>10);
            parser.addParameter('gaIterations',500,@(x)x>10);
            parser.addParameter('fminconIterations',500,@(x)x>10);
            parser.addParameter('verbose',true,@islogical);
            parser.addParameter('weighting',[],@isstruct);
            parser.parse(varargin{:});
            initial = parser.Results.initial;
            eps = parser.Results.eps;
            lb = parser.Results.lb;
            ub = parser.Results.ub;
            fix = parser.Results.fix;
            Usep = parser.Results.Usep;
            w = parser.Results.w;
            popSizeGA = parser.Results.gaPopulationSize;
            iterGA = parser.Results.gaIterations;
            iterFmincon = parser.Results.fminconIterations;
            verbose = parser.Results.verbose;
            weighting = parser.Results.weighting;
            if ~isstruct(eps)
                e = cell(1,length(msmrFieldnames));
                for k = 1:length(e); e{k} = eps; end
                eps = cell2struct(e,msmrFieldnames,2);
            end

            % Define optimization bounds. Ensure all parameters >0.
            for k = 1:length(msmrFieldnames)
                field = msmrFieldnames{k};
                if isfield(fix,field)
                    % No need to define upper and lower bounds for
                    % fixed (non-optimized) parameters.
                    continue;
                end
                if isfield(lb,field)
                    lb.(field) = ones(msmrFields.(field).len,1).*lb.(field);
                elseif ~isempty(initial)
                    lb.(field) = initial.(field).*(1-eps.(field));
                else
                    lb.(field) = zeros(msmrFields.(field).len,1);
                end
                if isfield(ub,field)
                    ub.(field) = ones(msmrFields.(field).len,1).*ub.(field);
                elseif ~isempty(initial)
                    ub.(field) = initial.(field).*(1+eps.(field));
                else
                    ub.(field) = Inf*ones(msmrFields.(field).len,1);
                end
                lb.(field) = max(0,lb.(field));
                ub.(field) = max(0,ub.(field));
            end
            
            % Ensure lithiation bounds remain between 0 and 1.
            for k = 1:length(msmrFieldnamesLith)
                field = msmrFieldnamesLith{k};
                if isfield(fix,field)
                    % No need to define upper and lower bounds for
                    % fixed (non-optimized) parameters.
                    continue;
                end
                lb.(field) = max(min(lb.(field),1),0);
                ub.(field) = max(min(ub.(field),1),0);
            end

            % Convert bounds to vectors.
            lb = ocp.MSMR.pack(lb,fix); ub = ocp.MSMR.pack(ub,fix);

            % Define the U1 > U2 > U3 > ... constraint. Ensures the
            % galleries do not "swap" positions so that constraints may be
            % applied to each individual gallery.
            orderU.Uj0 = toeplitz([-1 zeros(1,J-2)],[-1 1 zeros(1,J-2)])';
            fnames = setdiff(msmrFieldnames,{'Uj0'});
            for k = 1:length(fnames)
                field = fnames{k};
                orderU.(field) = zeros(msmrFields.(field).len,J-1);
            end
            A = ocp.MSMR.pack(orderU,fix)';
            B = -ones(J-1,1).*Usep(:);  % Enforce minimum gallery potential separation.

            % Define the sum(Xj)=1 constraint.
            sumX.Xj = ones(J,1);
            fnames = setdiff(msmrFieldnames,{'Xj'});
            for k = 1:length(fnames)
                field = fnames{k};
                sumX.(field) = zeros(msmrFields.(field).len,1);
            end
            Aeq = ocp.MSMR.pack(sumX,fix)';
            Beq = 1;

            % Specify partial initial population for the genetic algorithm.
            % (Currently only includes the initial guess if given).
            if ~isempty(initial)
                population = ocp.MSMR.pack(initial,fix)';
            else
                % No initial guess; let ga generate the initial population.
                population = [];
            end

            % Calculate the number of degrees of freedom in the optimization
            % (total number of MSMR variables to optimize).
            ndim = sum([ocp.MSMR.pack(msmrFields,fix).len]);

            % Configure hybrid Genetic Algorithm optimization.
            if verbose
                display = 'iter';
            else
                display = 'off';
            end
            optionsHybrid = optimoptions(@fmincon, ...
                'Display', display,...
                'MaxFunEvals', 1e6, ...
                'MaxIter', iterFmincon,...
                'TolFun', 1e-10, ...
                'TolX', 1e-10, ...
                'TolCon', 1e-10);
            options = optimoptions(@ga,...
                'Display', display, ...
                'PopulationSize', popSizeGA, ...
                'InitialPopulationMatrix', population, ...
                'MaxGenerations', iterGA, ...
                'FunctionTolerance', 1e-10,...
                'MaxStallGenerations', 100,...
                'HybridFcn', {@fmincon,optionsHybrid});

            [vect, Jcost] = ga(@cost,ndim,A,B,Aeq,Beq,lb,ub,[],options);
            est = ocp.MSMR.unpack(vect,fix);

            function Jcost = cost(vect)
                % Unpack parameter vector into MSMR model.
                model = ocp.MSMR.unpack(vect,fix);

                % Accumulate cost as differences between model prediction
                % and laboratory measurements (more than one
                % laboratory-derived OCP estimate may be supplied).
                Jcost = 0;
                for idxEstimate = 1:length(ocpestimates)
                    T = ocpestimates(idxEstimate).ocptest.temp;
                    vdata = ocpestimates(idxEstimate).V;
                    zdata = ocpestimates(idxEstimate).Z;
                    dzdata = max(zdata) - min(zdata);
                    dzdvData = ocpestimates(idxEstimate).dzdvRefV;
                    refVData = ocpestimates(idxEstimate).refV;
    
                    % We have to re-evalulate the model prediction for each
                    % OCP estimate, as the temperature T may differ.
                    [v, dzdv, z] = model.getOCP('volt',vmin,vmax,10000,T);
                    dzdv = abs(dzdv); % OCP estimates use absolute differential capacity!
    
                    % Convert model vectors over absolute lithiation to relative
                    % lithiation for comparison with data.
                    zmodel = (interp1(v,z,vdata,'linear','extrap') - model.zmin)/(model.zmax - model.zmin);
                    zmodel = min(zdata) + zmodel*dzdata;
                    dzdvModel = interp1(v,dzdv,refVData,'linear','extrap')*dzdata/(model.zmax - model.zmin);
    
                    % Compute residuals.
                    zResid = zmodel - zdata;
                    dzdvResid = dzdvModel - dzdvData;

                    % Apply weighting.
                    for weight = weighting(:)'
                        idx = weight.getInterval(vdata);
                        zResid(idx) = zResid(idx)*weight.multiplier;
                        idx = weight.getInterval(refVData);
                        dzdvResid(idx) = dzdvResid(idx)*weight.multiplier;
                    end
    
                    % Compute cost.
                    Jcost = Jcost + sum(zResid.^2) + w*sum(dzdvResid.^2);
                end
            end
        end

        function modAvg = average(mod1, mod2)
            %AVERAGE Average two MSMR models.

            if mod1.zmin ~= mod2.zmin || mod1.zmax ~= mod2.zmax
                error('Models have different absolute lithiation range');
            end

            % Parameters are already sorted by Uj0, so we average
            % corresponding parameters directly.
            modAvg = ocp.MSMR.unpack((ocp.MSMR.pack(mod1)+ocp.MSMR.pack(mod2))/2);
        end
    end
end