function lumped = lumpCellModel(std)
    %LUMPCELLMODEL Convert standard LMB model to lumped-parameter model.
    %
    % lumped = lumpCellModel(std) converts the pseudo-two-dimensional cell
    %   model structure, STD, to the lumped cell model structure, LUMPED.
    %
    % -- Changelog --
    % 2023.06.06 | Created | Wesley Hileman <whileman@uccs.edu>

    p = std; % ease notation in below computations

    % Convert P2D model to lumped form.
    l.const.F = p.const.F;
    l.const.R = p.const.R;
    l.const.T = p.const.T;
    l.const.Q = p.pos.epsilons*p.const.A*p.pos.L*p.const.F*p.pos.csmax*abs(p.pos.theta100-p.pos.theta0)/3600; % [Ah]
    l.const.psi = p.const.F*p.const.De*p.const.ce0/p.const.kappa/(1-p.const.t0plus)/p.const.T; % [V/K]
    l.const.kD = 2*p.const.R*(p.const.t0plus-1)*(1+p.const.dlnfdlnc)/p.const.F; % [V/K]
    l.const.Rc = p.const.Rc;
    l.pos.sigma = p.const.A*p.pos.sigma*(p.pos.epsilons)^p.const.brug/p.pos.L; % [1/Ω]
    l.pos.kappa = p.const.A*p.const.kappa*(p.pos.epsilone)^p.const.brug/p.pos.L; % [1/Ω]
    l.pos.qe = p.const.F*p.pos.epsilone*p.const.ce0*p.const.A*p.pos.L/(1-p.const.t0plus)/3600; % [Ah]
    l.pos.Dsref = p.pos.Dsref/p.pos.Rs^2; % [1/s]
    l.pos.nF = p.pos.nF; % [unitless]
    l.pos.nDL = p.pos.nDL; % [unitless]
    l.pos.Rf = p.pos.Rf/p.pos.as/p.const.A/p.pos.L; % [Ω]
    l.pos.Rdl = p.pos.Rdl/p.pos.as/p.const.A/p.pos.L; % [Ω]
    l.pos.Cdl = p.pos.Cdl*p.pos.as*p.const.A*p.pos.L; % [F]
    l.pos.U0 = p.pos.U0; % [V]
    l.pos.X = p.pos.X; % [unitless]
    l.pos.omega = p.pos.omega; % [unitless]
    l.pos.alpha = p.pos.alpha; % [unitless]
    l.pos.theta0 = p.pos.theta0;
    l.pos.theta100 = p.pos.theta100;
    l.pos.k0 = p.const.A*p.pos.L*p.pos.as*p.pos.k0; % [A]
    l.sep.kappa = p.const.A*p.const.kappa*(p.sep.epsilon)^p.const.brug/p.sep.L; % [1/Ω]
    l.sep.qe = p.const.F*p.sep.epsilon*p.const.ce0*p.const.A*p.sep.L/(1-p.const.t0plus)/3600; % [Ah]
    l.sep.nE = p.sep.nE;
    l.DL.kappa = p.const.A*p.const.kappa*(p.DL.epsilon)^p.const.brug/p.DL.L; % [1/Ω]
    l.DL.qe = p.const.F*p.DL.epsilon*p.const.ce0*p.const.A*p.DL.L/(1-p.const.t0plus)/3600; % [Ah]
    l.DL.nE = p.DL.nE;
    l.neg.k0 = p.neg.gamma*p.const.A*p.const.F*p.neg.k0*(p.neg.cs)^p.neg.alpha*p.const.ce0^(1-p.neg.alpha); % [A]
    l.neg.alpha = p.neg.alpha; % [unitless]
    l.neg.nDL = p.neg.nDL; % [unitless]
    l.neg.Rf = p.neg.Rf/p.neg.gamma/p.const.A; % [Ω]
    l.neg.Rdl = p.neg.Rdl/p.neg.gamma/p.const.A; % [Ω]
    l.neg.Cdl = p.neg.Cdl*p.neg.gamma*p.const.A; % [F]

    if isfield(l.pos,'U0')
        l.MSMR = true;
    else
        l.MSMR = false;
    end

    % Output assignment
    lumped = l;
end