function [R0, deltaV, phi_se0, phi_se2, phi_se3] = getPulseResistance(model,socPct,current,TdegC,method)
    %GETPULSERESISTANCE Compute pulse-resistance at a set of SOC and pulse-
    % current setpoints for a parameterized LMB cell. Fast implementation!
    %
    % GETPULSERESISTANCE(MODEL,SOC,CURRENT,TdegC) finds for the pulse-resistance of the
    %   LMB cell whose parameters are given by MODEL at the soc and
    %   current setpoints created by all pairs of values in vectors SOC and
    %   CURRENT. TdegC is the temperature is degC.
    %
    % GETPULSERESISTANCE(...,'exact') solves for pulse resistance using the exact
    %   solutions to the nonlinear equations. This is the default.
    %
    % GETPULSERESISTANCE(...,'linear') solves for pulse resistance using
    %   approximate linear solutions to the nonlinear equations. The
    %   resulting pulse-resistance does not vary with current magnitude.
    %   This may be useful for verification and rough approximations in
    %   which the overpotential remains small.
    %
    % Input:
    %   model   = PulseModel instance containing LMB parameters.
    %   socPct  = state-of-charge setpoints over which to evalulate R0 [%].
    %   current = pulse-current setpoints over which to evalulate R0 [A].
    %   method  = 'exact' or 'linear' (defaults to 'exact').
    %
    % Output: (Matricies where rows correspond to currents, columns to SOCs)
    %   R0      = matrix of pulse resistances [Ohm].
    %   deltaV  = matrix of cell voltage changes [V].
    %   phi_se0 = matrix of debiased phi_se(0,0+) values [V].
    %   phi_se2 = matrix of debiased phi_se(2,0+) values [V].
    %   phi_se3 = matrix of debiased phi_se(3,0+) values [V].
    %
    % Adapted from a version previously developed for LIB. With LMB, we 
    % no longer have an intercalation negative electrode, and need only 
    % solve an ODE in the positive electrode region (with bvp5c) and an 
    % algebraic equation at the surface of the negative electrode (with fzero).
    %
    % This function provides an initial guess to the solution of the
    % positive-electrode ODE by using the linear solution obtained by 
    % Taylor-series expansion of the reaction flux equation.
    %
    % Note: The built-in function `fzero` is extremely slow, so this
    % implementation uses the custom `fzeroFaster` function for improved 
    % performance.
    %
    % -- Changelog --
    % 01.06.2024 | Update for gen2 tookkit | Wesley Hileman
    % 09.10.2022 | Use linear approx to initialize bvp5c | Wesley Hileman
    % 08.16.2022 | Standardize model parameter names | Wesley Hileman
    % 03.18.2022 | Updated for LMB | 
    % Wesley Hileman <whileman@uccs.edu>
    % University of Colorado Colorado Springs

    % Initialize! ---------------------------------------------------------

    model = convertCellModel(model,'LLPM');
    soc = socPct/100;

    if nargin < 5
        method = 'exact';
    end
    if strcmp(method,'exact')
        useExact = true;
    elseif strcmp(method,'linear')
        useExact = false;
    else
        error('method input must be "exact" or "linear"');
    end

    % BVP solver options.
    % A little faster with "Vectorize" turned on.
    bvpOpts = bvpset('Stats','off','Vectorized','on');

    % Ensure we have column vectors.
    soc = soc(:);
    current = current(:);

    if isfield(model,'function')
        % Toolbox model.
        hasReg = @(reg)isfield(model.function,reg);
        getReg = @(reg)model.function.(reg);
        hasParam = @(reg,p)isfield(model.function,reg)&&isfield(model.function.(reg),p);
        getParam = @(reg,p)model.function.(reg).(p)(0,TdegC+273.15);
    else
        % Model parameters evalulated at setpoint.
        hasReg = @(reg)isfield(model,reg);
        getReg = @(reg)model.(reg);
        hasParam = @(reg,p)isfield(model,reg)&&isfield(model.(reg),p);
        getParam = @(reg,p)model.(reg).(p);
    end

    % Collect parameters.
    Rdlp = getParam('pos','Rdl');
    Rfp = getParam('pos','Rf'); 
    Xp = getParam('pos','X');
    U0p = getParam('pos','U0'); 
    omegap = getParam('pos','omega'); 
    alphap = getParam('pos','alpha'); 
    k0p = getParam('pos','k0'); 
    kappap = getParam('pos','kappa'); 
    sigmap = getParam('pos','sigma'); 
    theta0p = getParam('pos','theta0'); 
    theta100p = getParam('pos','theta100');
    Rdln = getParam('neg','Rdl'); 
    k0n = getParam('neg','k0'); 
    alphan = getParam('neg','alpha'); 
    if hasParam('const','Rc')
        % Rfn, 1/kappad, and 1/kappas appear in series in the model,
        % so we allow a single resistance parameter to be input instead
        % of three separate parameters.
        % Internally, we let 1/kappas "absorb" Rc and set kappad=Inf
        % (appears as 1/kappad=0 in the model equations) and Rfn=0.
        kappad = Inf;
        kappas = 1/getParam('const','Rc');
        Rfn = 0;
    else
        if hasReg('dll')
            kappad = getParam('dll','kappa'); 
        else
            % legacy model
            kappad = getParam('DL','kappa'); 
        end
        kappas = getParam('sep','kappa');
        Rfn = getParam('neg','Rf'); 
    end

    if Rfp == 0 && Rdlp == 0
        warning('Rdl and Rf (positive electrode) cannot both be zero');
        R0 = NaN(size(dv));
        return;
    end
    
    % Determine the OCP values corresponding to the SOC setpoints.
    electrode = MSMR(getReg('pos'));
    theta = theta0p + soc*(theta100p-theta0p);
    ocpData = electrode.ocp('theta',theta,'TdegC',TdegC,'npoints',100000);
    Uocp = ocpData.Uocp;
    
    % Precompute stuff to speed things up. Solvers evalulate equations 
    % many times, so anything we can precompute helps.
    f = electrode.f(TdegC);
    K1p = (1-alphap)*f;
    K2p = -alphap*f;

    % Allocate storage for phi_se values.
    phi_se0 = zeros([length(current), length(soc)]);
    phi_se2 = zeros(size(phi_se0));
    phi_se3 = zeros(size(phi_se0));

    % Solve! --------------------------------------------------------------
    
    % First, determine phi_se(0,0+) at the negative electrode for each 
    % applied current value. Note that phi_se(0,0+) is independent of SOC!
    Rctn = 1/k0n/f;
    for kc = 1:length(current)
        iapp = current(kc);
        guess = (Rctn*Rdln/(Rctn+Rdln) + Rfn)*iapp;  % Linearized solution!
        if useExact
            % Exact solution to nonlinear equation.
            phi_se0(kc,:) = fzeroFaster(@Neg_LiFluxBalanceEquation, guess);
        else
            % Linear approximation.
            phi_se0(kc,:) = guess;
        end
    end

    % Second, determine phi_se(2,0+) and phi_se(3,0+) in the positive
    % electrode for each SOC and applied current value.
    for kz = 1:length(Uocp)
        Urefp = Uocp(kz);
            
        % Compute fractions xj and i0 for each gallary (vectors)
        xp = Xp./(1+exp(f*(Urefp-U0p)./omegap)); % sum(xjp) = thetap
        i0p = k0p.*(xp.^(omegap.*alphap)).*((Xp-xp).^(omegap.*(1-alphap)));

        % Compute charge-transfer resistance for positive electrode
        % (we'll use it to guess an initial solution to the BVP below).
        Rctp = 1/f/sum(i0p);

        % Linearized resistance of the interface for use in guessing
        % solutions to the BVP.
        Rinterfacep = Rfp+(Rctp*Rdlp)/(Rctp+Rdlp);
        
        % For use in solving the reaction flux equation for If+dl inside
        % the ODE function defined below.
        posFn = @(Ifdl,phise)(...
            sum(i0p.*(exp(K1p.*(phise-Rfp*Ifdl))-exp(K2p.*(phise-Rfp*Ifdl)))) ...  % If
            + (phise-Rfp*Ifdl)/Rdlp ... % Idl
            - Ifdl ...
        );
        
        for kc = 1:length(current)
            ik = current(kc);

            % Solve the ODE! ----------------------------------------------

            % Initial guess at the solution. This is the result when we
            % linearize the reaction flux equation w/r/t the overpotential
            % by Taylor series expansion. Good approximation when the
            % overpotential is small.
            zeta = sqrt((1/sigmap+1/kappap)/Rinterfacep);
            %A = [
            %    zeta            -zeta
            %    zeta*exp(zeta)  -zeta*exp(-zeta)
            %];
            % Use 2x2 inverse formula directly (better numerical
            % properties as zeta->inf).
            Ainv = [
                -exp(-zeta)     1
                -exp(zeta)      1
            ]/zeta/(exp(zeta)-exp(-zeta));
            b = [
                +ik/kappap
                -ik/sigmap
            ];
            coeff = Ainv*b;
            phise2 = coeff(1)+coeff(2);
            phise3 = coeff(1)*exp(zeta)+coeff(2)*exp(-zeta);
            if abs(phise2)>abs(phise3)
                phiseMax = phise2;
            else
                phiseMax = phise3;
            end
            ifdlMax = phiseMax/Rinterfacep;
            etaMax = phiseMax - ifdlMax*Rfp;
            yinit = @(x) [
                coeff(1)*exp(zeta*x) + coeff(2)*exp(-zeta*x)
                coeff(1)*zeta*exp(zeta*x) - coeff(2)*zeta*exp(-zeta*x)
            ]; 

            % NOTE: Okay to use linear solution if overpotential is small.
            % This speeds up the computation time in some cases.
            if useExact && 0.5*f*abs(etaMax)>0.01
                % Find the exact solution to the ODE.
                init.solver = 'bvpinit';
                init.x      = [0 0.5 1];
                init.y      = yinit([0 0.5 1]);
                init.yinit  = yinit;
                sol5c = bvp5c(@odefun,@bcfun,init,bvpOpts);
                phi_se2(kc,kz) = sol5c.y(1,1);
                phi_se3(kc,kz) = sol5c.y(1,end);
            else
                % Use the approximate solution to the ODE.
                y0 = yinit(0); y1 = yinit(1);
                phi_se2(kc,kz) = y0(1);
                phi_se3(kc,kz) = y1(1);
            end
        end
    end

    % Third, compute deltaV and the pulse resistance.
    Iapp = repmat(current,1,length(soc));
    deltaV = -(1/kappad + 1/kappas + 1/(sigmap+kappap))*Iapp ...
             - phi_se0 ...
             + sigmap/(kappap+sigmap)*phi_se2 ...
             + kappap/(kappap+sigmap)*phi_se3;
    R0 = -deltaV./Iapp;

    % Solver Functions ----------------------------------------------------
    % NOTE: these use some of the variables defined above.

    function residual = Neg_LiFluxBalanceEquation(phi_se)
        %NEG_LIFLUXBALANCEEQUATION For use in solving the nonlinear 
        % negative-electrode flux density equation for phi_se with fzero.

        eta = phi_se - Rfn*iapp;
        If = k0n*(exp((1-alphan)*f*eta) - exp(-alphan*f*eta));
        Ifdl = If + eta/Rdln;
        residual = iapp - Ifdl;
    end

    % Define ODE function: "x" (spatial position) input is not used
    % y    = [phisep;phisep']
    % dydx = [phisep';phisep'']
    function dydx = odefun(~,y)
        dydx = zeros(size(y)); 
        
        % First we get d(phi_se)/dx from the second element in y.
        dydx(1,:) = y(2,:);
        
        % Now we solve for d^2(phi_se)/dx^2 using phi_se, which is the
        % first element of y.
        if Rfp == 0 && Rdlp ~= 0
            dydx(2,:) = (1/sigmap+1/kappap)...
                *(i0p*(exp(K1p*y(3,:))-exp(K2p*y(3,:)))...
                +y(1,:)/Rdlp);
        elseif Rfp ~= 0 && Rdlp == 0
            dydx(2,:) = (1/sigmap+1/kappap)*y(1,:)/Rfp;
        else
            % 08.18.20222 | WH
            % Better initial guess reduces fmincon failures;
            % choosing ifdl=phi/Rfp makes eta=0, means maximum residual
            % is ifdl=phi/Rfp itself.
            for k = 1:size(y,2)
                phi = y(1,k);
                if zeta < 1
                    guessP = -ik;
                else
                    guessP = phi/Rfp;
                end
                try
                    ifdl = fzeroFaster(posFn,guessP,[],phi);
                catch
                    error('could not find if+dl (1)')
                end
                if isnan(ifdl)
                    error('could not find if+dl (2)');
                end
                dydx(2,k) = (1/sigmap+1/kappap)*ifdl;
            end
        end
    end

    % Set up boundary conditions
    % y = [phisen;phisen';phisep;phisep'];
    % yl is the left-side boundary, yr is the right-side boundary
    function res = bcfun(yl,yr)
        res = [yl(2)-ik/kappap;  % xtilde = 2
               yr(2)+ik/sigmap]; % xtilde = 3
    end
end