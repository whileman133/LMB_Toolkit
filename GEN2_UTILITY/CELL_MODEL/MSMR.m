classdef MSMR < handle
    %MSMR Multiple-site multiple-reaction (MSMR) OCP model for a porous
    %  electrode.
    %
    % -- Usage --
    %
    % electrode = MSMR(params,'pos') constructs an MSMR model from the
    %   parameter structure PARAMS for a positive porous electrode
    %   (thetamin=theta100 and thetamax=theta0). PARAMS is a structure of
    %   parameters describing the positive electrode. The parameter values
    %   may either be functions (in which case they will be evalulated at
    %   the default SOC and temperature, which can be changed by keyword
    %   arguments) or numeric values. The following vector parameters are 
    %   required: U0, X, and omega. The following scalar parameters are 
    %   also required: theta0 and theta100.
    %
    % electrode = MSMR(params,'neg') constructs an MSMR model from the
    %   parameter structure PARAMS for a negative porous electrode
    %   (thetamin=theta0 and thetamax=theta100).
    %
    % electrode = MSMR(...,'TdegC',TdegC) sets the temperature setpoint to 
    %   TdegC instead of the default 25degC. This option is used only to
    %   evalulate function-valued parameters.
    %
    % electrode = MSMR(...,'name',name) assigns a name to the porous 
    %   electrode model.
    %
    % electrode = MSMR(...,'reference',reference) assigns a literature 
    %   reference to the porous electrode model.
    %
    % -- Built-in models --
    %
    % electrode = MSMR.C6() constructs an MSMR model for lithiated graphite. 
    %   Mark Verbrugge et al. 2017 J. Electrochem. Soc. 164 E3243.
    %
    % electrode = MSMR.NMC622() constructs an MSMR model for NMC622. 
    %   Mark Verbrugge et al. 2017 J. Electrochem. Soc. 164 E3243.
    %
    % electrode = MSMR.LMO() constructs an MSMR model for LMO.
    %   Daniel R. Baker and Mark W. Verbrugge 2021 J. Electrochem. Soc. 168 050526.
    %
    % -- Methods --
    %
    %  electrode.ocp() compute OCP versus lithiation.
    %  electrode.Rct() compute charge-transfer resistance versus lithiation.
    %  electrode.lith2ocp() compute OCP at a specific lithiation point.
    %  electrode.plotOCP() plot the OCP, OCP slope, and OCP curvature for
    %      quick inspection.
    %
    %  Type help MSMR.<method> for additional information on a method.
    %
    % -- Changelog --
    % 07.02.2023 | Single output structs for ocp(), Rct() | Wesley Hileman
    % 06.23.2023 | Depricate 'electrode' parameter | Wesley Hileman
    % 01.23.2023 | Simplify usage | Wesley Hileman
    % 02.05.2022 | Created | Wesley Hileman <whileman@uccs.edu>

    properties(Constant)
        R = 8.3144598;      % Molar gas constant [J/mol K]
        F = 96485.3329;     % Faraday constant [C/mol]
    end

    properties(SetAccess=protected)
        J        % Number of galleries.
        Xj       % Fraction of sites belonging to each gallery.
        Uj0      % Standard electrode potentials for each reaction.
        Wj       % Nernst equation shape factor for each gallery.
        zmin     % Minimum lithiation.
        zmax     % Maximum lithiation.
        name     % Model name (if applicable).
        ref      % Citation/reference (if applicable).
        sortidx  % Set when 'sortParams' option is true for keeping track
                 % of sort order of MSMR parameters. 
    end

    methods
        function obj = MSMR(params, varargin)
            %MSMR Construct an MSMR object with the specified parameters.
            %  Xj, Uj0, Wj, zmin, zmax

            p = inputParser;
            p.addRequired('params',@(x)isstruct(x));
            p.addParameter('TdegC',25);
            p.addParameter('name','');
            p.addParameter('reference','');
            p.addParameter('sortParams',false,@islogical);
            p.parse(params, varargin{:});
            TdegC = p.Results.TdegC;
            T = TdegC+273.15;
            name = p.Results.name;
            reference = p.Results.reference;
            sortParams = p.Results.sortParams;

            if all(isfield(params,{'theta0','theta100'}))
                zmin = params.theta100;
                zmax = params.theta0;
            elseif all(isfield(params,{'thetamin','thetamax'}))
                zmin = params.thetamin;
                zmax = params.thetamax;
            elseif all(isfield(params,{'zmin','zmax'}))
                zmin = params.zmin;
                zmax = params.zmax;
            else
                error(['Could not identify thetamin and thetamax from ' ...
                    'supplied parameters!']);
            end
            if isa(zmin,'function_handle'), zmin = zmin(0,T); end
            if isa(zmax,'function_handle'), zmax = zmax(0,T); end
            if zmin>zmax
                tmp = zmin;
                zmin = zmax;
                zmax = tmp;
            end

            Xj = params.X;
            Uj0 = params.U0;
            Wj = params.omega;
            if isa(Xj,'function_handle'), Xj = Xj(0,T); end
            if isa(Uj0,'function_handle'), Uj0 = Uj0(0,T); end
            if isa(Wj,'function_handle'), Wj = Wj(0,T); end

            obj.Xj = Xj(:);
            obj.Uj0 = Uj0(:);
            obj.Wj = Wj(:);
            obj.J = length(Xj);
            obj.zmin = zmin;
            obj.zmax = zmax;
            obj.name = name;
            obj.ref = reference;

            % Sort parameters by equilibrium potential.
            obj.sortidx = [];
            if sortParams
                [obj.Uj0, idx] = sort(obj.Uj0,'ascend');
                obj.Xj = obj.Xj(idx);
                obj.Wj = obj.Wj(idx);
                obj.sortidx = idx;
            end
        end

        function ocpData = ocp(obj, varargin)
            %OCP Calculate OCP predicted by the MSMR model.
            %
            % ocpData = OCP() computes the MSMR OCP curve 
            %   over the default lithiation range (thetamin to thetamax).
            %   Potential points are evenly spaced.
            %
            % [...] = OCP('vmin',vmin,'vmax',vmax) computes the MSMR
            %   OCP curve over the potential range VMIN to VMAX.
            %   Potential points are evenly spaced.
            %
            % [...] = OCP('thetamin',tmin,'thetamax',tmax) computes 
            %   the MSMR OCP curve over the lithiation range TMIN to TMAX.
            %   Potential points are evenly spaced. Lithiation points are
            %   *unevely* spaced.
            %
            % [...] = OCP('voltage',voltage) computes the MSMR
            %   OCP curve at potential points given by the vector VOLTAGE.
            %
            % [...] = OCP('theta',theta) computes the MSMR
            %   OCP curve at lithiation points given by the vector THETA.
            %   NPOINTS determines the accuracy of the approximation.
            %
            % [...] = OCP(...,'TdegC',T) computes the OCP at temperature 
            %   T degC.
            %
            % [...] = OCP(...,'npoints',npoints) uses NPOINTS points to
            %   generate the potential vector (if the potential vector is 
            %   not explicitly supplied, but rather specified by min/max 
            %   bounds in either potential or lithiation). The default
            %   number of points is 10000.
            %
            % Output structure:
            %   ocpData.Uocp      OCP vector
            %   ocpData.dUocp     OCP slope vector
            %   ocpData.d2Uocp    OCP curvature vector
            %   ocpData.theta     Lithiation vector
            %   ocpData.xj        Lithiation matix: d1=gallery d2=SOC-point
            %   ocpData.U0        MSMR U0 vector
            %   ocpData.X         MSMR X vector
            %   ocpData.omega     MSMR omega vector
            %   ocpData.thetamin  Min. lithiation of the electrode
            %   ocpData.thetamax  Max. lithiation of electrode

            p = inputParser;
            p.addOptional('vmin',[]);
            p.addOptional('vmax',[]);
            p.addOptional('voltage',[]);
            p.addOptional('thetamin',[]);
            p.addOptional('thetamax',[]);
            p.addOptional('theta',[]);
            p.addOptional('TdegC',25);
            p.addOptional('npoints',10000);
            p.parse(varargin{:});
            vmin = p.Results.vmin;
            vmax = p.Results.vmax;
            voltage = p.Results.voltage;
            thetamin = p.Results.thetamin;
            thetamax = p.Results.thetamax;
            theta = p.Results.theta;
            TdegC = p.Results.TdegC;
            npoints = p.Results.npoints;

            % Ensure everything is a row vector.
            voltage = voltage(:)';
            theta = theta(:)';

            % Determine potential vector.
            if ~isempty(vmin) && ~isempty(vmax)
                Uocp = linspace(vmin, vmax, npoints);
            elseif ~isempty(voltage)
                Uocp = voltage;
            elseif ~isempty(thetamin) && ~isempty(thetamax)
                vmax = obj.lith2ocp(thetamin,'TdegC',TdegC);
                vmin = obj.lith2ocp(thetamax,'TdegC',TdegC);
                Uocp = linspace(vmin, vmax, npoints);
            elseif ~isempty(theta)
                vmax = obj.lith2ocp(min(theta),'TdegC',TdegC);
                vmin = obj.lith2ocp(max(theta),'TdegC',TdegC);
                Uocp = linspace(vmin, vmax, npoints);
            else
                vmax = obj.lith2ocp(obj.zmin,'TdegC',TdegC);
                vmin = obj.lith2ocp(obj.zmax,'TdegC',TdegC);
                Uocp = linspace(vmin, vmax, npoints);
            end

            % Compute stoichiometry vector.
            f = obj.f(TdegC);
            gj = exp(f*(Uocp-obj.Uj0)./obj.Wj);
            xj = obj.Xj./(1+gj);
            Z = sum(xj,1);

            % Compute differential OCP.
            dUocp = -1./(f*sum((obj.Xj./obj.Wj).*gj./(1+gj).^2,1));
            d2Uocp = (dUocp).^3.*...
                     sum((obj.Xj./obj.Wj.^2).*gj.*(1-gj)./(1+gj).^3,1)...
                     *f^2;

            % Interpolate over requested stiochiometry vector (if supplied).
            if ~isempty(theta)
                Unew = interp1(Z,Uocp,theta,'linear','extrap');
                dUdzNew = interp1(Z,dUocp,theta,'linear','extrap');
                d2UocpNew = interp1(Z,d2Uocp,theta,'linear','extrap');
                xjNew = zeros(obj.J,length(Unew));
                for j = 1:obj.J
                    xjNew(j,:) = interp1(Uocp,xj(j,:),Unew,'linear','extrap');
                end
                Uocp = Unew;
                dUocp = dUdzNew;
                d2Uocp = d2UocpNew;
                xj = xjNew;
            else
                theta = Z;
            end

            % Collect output data.
            ocpData.Uocp = Uocp;
            ocpData.dUocp = dUocp;
            ocpData.d2Uocp = d2Uocp;
            ocpData.theta = theta;
            ocpData.xj = xj;
            ocpData.TdegC = TdegC;
            ocpData.f = f;
            ocpData.npoints = npoints;
            ocpData.U0 = obj.Uj0;
            ocpData.omega = obj.Wj;
            ocpData.X = obj.Xj;
            ocpData.thetamin = obj.zmin;
            ocpData.thetamax = obj.zmax;
            ocpData.origin__ = 'MSMR.ocp';
        end

        function ctData = Rct(obj, params, varargin)
            %RCT Compute the MSMR charge-transfer resistance versus SOC.
            %
            % -- Usage --
            % ctData = RCT(params,...) computes the MSMR charge-transfer
            %   resistance at the SOC setpoints specified in (...) using
            %   exactly the same input syntax as MSMR.OCP(). PARAMS is
            %   a structure containing vector fields specifying the 
            %   kinetics parameters k0 and alpha.
            %
            % Output structure: Everying returned by MSMR.ocp plus:
            %   ctData.Rct   Charge-transfer (CT) resistance vector
            %   ctData.i0    Exchange-current (EC) vector
            %   ctData.Rctj  CT resistance matrix: d1=gallery d2=SOC-point
            %   ctData.i0j   EC matrix: d1=gallery d2=SOC-point
            %
            % This function assumes theta_e=1 to perform the computation.
            % (Holds for small-signal linear approximation used in transfer 
            % function model).
            
            ocpData = obj.ocp(varargin{:});
            ctData = RctCachedOCP(obj, params, ocpData);
        end

        function ctData = RctCachedOCP(obj, params, ocpData)
            % RCTCACHEDOCP Same as Rct(), but uses provided ocpData instead
            % of calling obj.OCP(); allows for code optimization.

            xj = ocpData.xj;
            TdegC = ocpData.TdegC;
            T = TdegC+273.15;

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

            % Make sure to sort k0, alpha in same order as U, X, omega!
            if ~isempty(obj.sortidx)
                k0 = k0(obj.sortidx);
                alpha = alpha(obj.sortidx);
            end

            % Now compute the charge-transfer resistance.
            i0j = k0.*xj.^(omega.*alpha).*(X-xj).^(omega.*(1-alpha));
            i0 = sum(i0j,1);
            Rct = 1./i0/ocpData.f;

            % Collect output data.
            ctData = ocpData;
            ctData.Rct = Rct;
            ctData.i0 = i0;
            ctData.Rctj = 1./i0j./ocpData.f;
            ctData.i0j = i0j;
            ctData.origin__ = 'MSMR.Rct';
        end

        function [Uocp, dUocp, xj, data] = lith2ocp(obj, z, varargin)
            %LITH2OCP Calculate the OCP corresponding to a
            %  particular surface lithiation by solving the nonlinear
            %  equation in Uocp. This implementation uses fzero.
            %
            % Input:
            %   z     = surface lithiation / stiochiometry
            %   TdegC = temperature, defaults to 25degC. Supply as a
            %           keyword ('name',value) optional argument.
            %
            % Output:
            %   Uocp   = open-circuit potential corresponding to z
            %   dUocp  = d(Uocp)/d(theta) at the specified lithiation
            %   xj     = partial surface stiochiometries corresponding to z
            %            and Uocp

            p = inputParser;
            p.addOptional('TdegC',25);
            p.parse(varargin{:});
            TdegC = p.Results.TdegC;
            f = obj.f(TdegC);

            Uocp = fzero(@lithiationBalanceEquation, mean(obj.Uj0));
            gj = exp(f*(Uocp-obj.Uj0)./obj.Wj);
            xj = obj.Xj./(1+gj);
            dxjdU = -f.*obj.Xj.*gj./obj.Wj./(1+gj).^2;
            dUocp = 1./sum(dxjdU,1);

            data.gj = gj;
            
            function residual = lithiationBalanceEquation(U)
                residual = z - sum(obj.Xj./(1+exp(f*(U-obj.Uj0)./obj.Wj)));
            end
        end

        function figs = plotOCP(obj, varargin)
            data = obj.ocp(varargin{:});
            Uocp = data.Uocp;
            dUocp = data.dUocp;
            theta = data.theta;

            % Plot formatter function.
            function format
                if exist('thesisFormat','file')
                    thesisFormat;
                end
            end

            colors = [cool(obj.J); 0 0 0];
            galleryLabels = arrayfun( ...
                @(j)sprintf('j = %d',j),1:obj.J, ...
                'UniformOutput',false);

            figs.ocp = figure; 
            colororder(colors);
            plot(data.xj,Uocp,':'); hold on;
            plot(theta,Uocp);
            xlabel('Fractional Lithiation, \theta_s');
            ylabel('Potential vs. Li/Li^+, U_{ocp} [V]');
            title(sprintf('%s: OCP',obj.name));
            legend(galleryLabels{:},'\Sigma');
            format;

            figs.docp = figure;
            plot(theta,dUocp);
            xlabel('Fractional Lithiation, \theta_s');
            ylabel('dU_{ocp}/d\theta [V]');
            title(sprintf('%s: OCP Slope',obj.name));
            format;

            figs.d2ocp = figure;
            plot(theta,data.d2Uocp);
            xlabel('Fractional Lithiation, \theta_s');
            ylabel('d^2U_{ocp}/d\theta^2 [V]');
            title(sprintf('%s: OCP Curvature',obj.name));
            format;
        end
    end

    methods(Static)
        function f = f(T)
            %F Compute f=F/(R*T) for the supplied T in [C].

            if ~exist('T','var')
                T = 25;
            end
            F = MSMR.F;
            R = MSMR.R;
            f = F/R/(T+273.15);
        end

        function obj = C6()
            %C6 construct a model for lithiated graphite from Verburgge
            %  2017.
            params.X = [0.43336 0.23963 0.15018 0.05462 0.06744 0.05476];
            params.U0 = [0.08843 0.12799 0.14331 0.16984 0.21446 0.36325];
            params.omega = [0.08611 0.08009 0.72469 2.53277 0.09470 5.97354];
            params.theta0 = 0.01;
            params.theta100 = 0.99;
            obj = MSMR(params, ...
                'name','C6 Electrode', ...
                'reference',['Mark Verbrugge et al. 2017 J. ' ...
                             'Electrochem. Soc. 164 E3243'], ...
                'sortParams',true);
        end

        function obj = NMC622()
            %NMC622 Construct a model for NMC622 from Verbrugge 2017.
            params.X = [0.13442 0.3246 0.21118 0.32980];
            params.U0 = [3.62274 3.72645 3.90575 4.22955];
            params.omega = [0.96710 1.39712 3.505 5.52757];
            params.theta0 = 0.99;
            params.theta100 = 0.1325;
            obj = MSMR(params, ...
                'name','NMC622 Electrode', ...
                'reference',['Mark Verbrugge et al. 2017 J. ' ...
                             'Electrochem. Soc. 164 E3243'], ...
                'sortParams',true);
        end

        function obj = LMO()
            %LMO Construct a model for LMO from Verbrugge 2021
            params.X = [1-0.60331 0.60331];
            params.U0 = [4.16756 4.02477];
            params.omega = [1.12446 1.71031];
            params.theta0 = 0.99;
            params.theta100 = 0.01;
            obj = MSMR(params, ...
                'name','LMO Electrode', ...
                'reference',['Daniel R. Baker and Mark W. Verbrugge 2021 ' ...
                             'J. Electrochem. Soc. 168 050526'], ...
                'sortParams',true);
        end
    end
end
