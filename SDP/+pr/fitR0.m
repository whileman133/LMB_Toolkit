function [result, modelspec, trajectory] = fitR0(surface,msmr,lb,ub,varargin)
    %FITR0 Regress a pulse-resistance model to a resistance surface.

    % Define parameters of the pulse model.
    params.pos.U0        = fastopt.param('fix',msmr.Uj0);
    params.pos.X         = fastopt.param('fix',msmr.Xj);
    params.pos.omega     = fastopt.param('fix',msmr.Wj);
    params.pos.theta0    = fastopt.param('fix',msmr.zmax);
    params.pos.theta100  = fastopt.param('fix',msmr.zmin);
    params.pos.alpha     = fastopt.param('fix',0.5*ones(msmr.J,1));
    params.pos.k0        = fastopt.param('len',msmr.J,'logscale',true);
    params.pos.Rf        = fastopt.param;
    params.pos.Rdl       = fastopt.param;
    params.pos.kappa     = fastopt.param;
    params.pos.sigma     = fastopt.param;
    params.neg.alpha     = fastopt.param('fix',0.5);
    params.neg.k0        = fastopt.param('logscale',true);
    params.neg.Rdl       = fastopt.param;
    % kappa(d), kappa(s), and Rf(n) are not separately identifiable, so
    % we lump them into the cell contact/tab resistance Rc.
    params.const.Rc      = fastopt.param;
    modelspec = fastopt.modelspec(params);

    % Convert lower/upper bounds to vectors.
    lb = fastopt.pack(lb,modelspec,'coerce',true); 
    ub = fastopt.pack(ub,modelspec,'coerce',true);

    % Extract pulse-resistance from surface object.
    iapp = surface.iapp;
    soc = surface.z;
    QAh = surface.QAh;
    TdegC = surface.temp;
    if isfield(surface,'R0')
        R0 = surface.R0;
        R0bound = surface.R0bound;
    else
        [R0, R0bound] = surface.getParam(@(ecm)ecm.R0);
    end

    % Throw out all pulse-resistance estimates at rates between +1C and -1C.
    % (lab data has poor SNR in this range).
    idx1C = (-QAh<=iapp)&(iapp<=+QAh);
    iapp(idx1C) = [];
    R0(idx1C,:) = [];
    R0bound(idx1C,:) = [];

    % Perform the regression.
    [result, trajectory] = ...
        fastopt.particleswarm(@cost,modelspec,lb,ub,varargin{:});
    
    function J = cost(model)
        R0model = pr.getPulseResistance(model,soc,iapp,TdegC,'exact');
        residuals = R0model - R0;
    
        % Compute cost as sum of square residuals weighted by 3sigma bounds.
        weightedResiduals = residuals./R0bound;
        J = sum(weightedResiduals(:).^2);
    end
end

