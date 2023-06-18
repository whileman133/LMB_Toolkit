% function [zk,boundzk,ekfData,Xind] = iterSPKF(vk,ik,Tk,spkfData)
% 
% Inputs:
%   vk       = The measured cell voltage for this time sample (V)
%   ik       = The measured cell current for this time sample (A)
%   Tk       = The measured cell temperature for this time sample (degC)
%   spkfData = Data structure initialized with initKF and updated here
% Outputs:
%   zk       = Comprises three sections, stacked vertically. The top of zk
%              holds all output variables (with nonlinear corrections
%              applied) in the same order they are organized in the ROMs.
%              Under this, zk has a (scalar) voltage estimate. Under this,
%              zk has a (scalar) SOC estimate (0..1).
%   boundzk  = A vector of 3-sigma confidence bounds for every element in
%              zk. Same dimension as zk.
%   spkfData = Data structure used by SPKF, with updated contents.
%   Xind     = Optional output, containing indices and scaling factors for
%              the individual ROMs used in the blend (this should not be
%              needed by most users so is a utility/debug output only).
%
% This function updates the SPKF based on measurements made during one
% sampling period. It outputs estimates of internal variables, voltage, and
% SOC, as well as 3-sigma confidence bounds on those estimates.
%
% Copyright (c) 2021 by Gregory L. Plett and M. Scott Trimboli of the
% University of Colorado Colorado Springs (UCCS). This work is licensed
% under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 Intl.
% License, v. 1.0. It is provided "as is", without express or implied
% warranty, for educational and informational purposes only.
%
% This file is provided as a supplement to: Plett, Gregory L. and Trimboli,
% M. Scott, "Battery Management Systems, Volume III, Physics-Based
% Methods," Artech House, 2021. 

function [zk,boundzk,spkfData,Xind] = iterSPKF(vk,ik,Tk,spkfData)
  persistent ind loc numT numZ V Z Tpts Zpts n warnCount cellData F R 
  if isempty(ind)
    checkBlendMethod;
    setupIndsLocs;
    % All temperature and SOC setpoints (sort in ascending order)
    Tpts = unique(spkfData.T);
    Zpts = unique(spkfData.Z);
    numT = length(Tpts); numZ = length(Zpts);
    switch spkfData.method
      case 'OB'
        V = zeros(1,8*spkfData.n+3); % voltage sigma points
        Z = zeros(spkfData.nz+2,8*spkfData.n+3); % output-variable sigma points
      case 'MB'
        V = zeros(1,2*spkfData.n+3); % voltage sigma points
        Z = zeros(spkfData.nz+2,2*spkfData.n+3); % output-variable sigma points
    end
    cellData  = spkfData.cellData;
    F         = cellData.const.F;
    R         = cellData.const.R;     
    n = spkfData.n;
    spkfData.ind = ind; % helpful to save for later when plotting
    spkfData.loc = loc;
    spkfData.Q = cellData.function.const.Q(); % cell capacity in Ah    
    warnCount = 0;
    
    if isempty(spkfData.SOC0)
      ocvFn = @(x) vk - (cellData.function.pos.Uocp(cellData.function.pos.soc(x)) ...
                 - cellData.function.neg.Uocp(cellData.function.neg.soc(x)));
      spkfData.SOC0 = fzero(ocvFn,0.5); % guess at initial SOC by voltage 
      fprintf('Initializing SOC0 to %g%%\n',spkfData.SOC0*100);
    end    
  end

  if warnCount > spkfData.maxWarn % SPKF is probably broken/lost
    zk = NaN(spkfData.nz+2,1); boundzk = zk;    
    Xind.gamma = NaN(4,1); Xind.theT = NaN(4,1); Xind.theZ = NaN(4,1);
    return
  end
  
  % Convert degC to K
  if Tk > 100
    warning('iterSPKF assumes that Tk is in Celsius. Converting...');
  else
    Tk = Tk + 273.15; % Convert to Kelvin...
  end
  
  % ----------------------------------------------------------------------
  % Steps 1a and 1b... update state(s) prediction and coviariance
  switch spkfData.method
    case 'OB'
      % Need to update state and covariance for ALL models
      for theT = 1:numT
        for theZ = 1:numZ
          A = spkfData.M(theT,theZ).A;
          xhat = spkfData.M(theT,theZ).xhat;
          xhat = A.*xhat + spkfData.priorI; % don't forget to use prior current!
          spkfData.M(theT,theZ).xhat = xhat;

          SigmaX = spkfData.M(theT,theZ).SigmaX;
          SigmaX = diag(A)*SigmaX*diag(A) + spkfData.SigmaW;
          spkfData.M(theT,theZ).SigmaX = SigmaX;
        end
      end
      spkfData.x0      = spkfData.x0 + spkfData.priorI; % not ik!!
      spkfData.SigmaX0 = spkfData.SigmaX0 + spkfData.SigmaW;

      % Compute time-update prediction of SOC
      SOC = spkfData.SOC0 - spkfData.x0*(spkfData.Ts/(3600*spkfData.Q)); 
    case 'MB'
      % Need to update state and covariance of single model
      SOC = spkfData.SOC0 - spkfData.xhat(end)*(spkfData.Ts/(3600*spkfData.Q)); 
      Xind = getXind(Tk,SOC);
      A1 = spkfData.M(Xind.theT(1),Xind.theZ(1)).A;
      A2 = spkfData.M(Xind.theT(2),Xind.theZ(2)).A;
      A3 = spkfData.M(Xind.theT(3),Xind.theZ(3)).A;
      A4 = spkfData.M(Xind.theT(4),Xind.theZ(4)).A;
      AMB = [[A1,A2,A3,A4]*Xind.gamma;1]; % model-blend A
      spkfData.xhat = AMB.*spkfData.xhat + spkfData.priorI;
      
      spkfData.SigmaX = diag(AMB)*spkfData.SigmaX*diag(AMB) + spkfData.SigmaW;
      SOC = spkfData.SOC0 - spkfData.xhat(end)*(spkfData.Ts/(3600*spkfData.Q)); 
  end

  % ----------------------------------------------------------------------
  % Step 1c... predict measured voltage
  % output blend: operates only on four nearest neighbors
  %               height of X and Xtilde are both = 4*n+1; width is
  %               number of sigma pts 
  % model blend:  height of X and Xtilde are both n+1; width is number
  %               of sigma points
  Xind = getXind(Tk,SOC);
  [X,Xtilde] = getSigmaPoints(Xind);  
  for theX = 1:size(X,2)
    V(theX) = getVariables(X(:,theX),ik,Xind,Tk); % this uses ik, not priorI
  end
  vhat = V*spkfData.alpham;
  
  % ----------------------------------------------------------------------
  % Step 2a... operates only on four nearest neighbors
  Vtilde = vhat - V;
  SigmaVtilde = Vtilde * diag(spkfData.alphac) * Vtilde' + spkfData.SigmaV;
  SigmaXVtilde = Xtilde * diag(spkfData.alphac) * Vtilde';
  L = SigmaXVtilde/SigmaVtilde;
  
  % ----------------------------------------------------------------------
  % Steps 2b and 2c... operate only on four nearest neighbors
  residual = vk - vhat;
  switch spkfData.method
    case 'OB'
      for theModel = 1:4
        xhat = spkfData.M(Xind.theT(theModel),Xind.theZ(theModel)).xhat;
        Lind = n*(theModel-1)+1:n*theModel;
        xhat = xhat + L(Lind)*residual;
        spkfData.M(Xind.theT(theModel),Xind.theZ(theModel)).xhat = xhat;

        SigmaX = spkfData.M(Xind.theT(theModel),Xind.theZ(theModel)).SigmaX;
        SigmaX = SigmaX - L(Lind)*SigmaVtilde*L(Lind)';
        [~,SS,VV] = svd(SigmaX);
        HH = VV*SS*VV';
        SigmaX = (SigmaX + SigmaX' + HH + HH')/4; % Help maintain robustness

        % Q-bump code
        if residual^2>9*SigmaVtilde % bad voltage estimate by 3-SigmaX, bump Q 
          fprintf('Increasing SigmaX\n');
          SigmaX = SigmaX*2;
        end

        spkfData.M(Xind.theT(theModel),Xind.theZ(theModel)).SigmaX = SigmaX;
      end
      spkfData.x0 = spkfData.x0 + L(end)*residual;
      spkfData.SigmaX0 = spkfData.SigmaX0 - L(end)*SigmaVtilde*L(end);
      
      % Compute measurement-update estimate of SOC
      SOC = spkfData.SOC0 - spkfData.x0*(spkfData.Ts/(3600*spkfData.Q));         
    case 'MB'
      spkfData.xhat = spkfData.xhat + L*residual;
      SigmaX = spkfData.SigmaX;
      SigmaX = SigmaX - L*SigmaVtilde*L';
      [~,SS,VV] = svd(SigmaX);
      HH = VV*SS*VV';
      SigmaX = (SigmaX + SigmaX' + HH + HH')/4; % Help maintain robustness

      % Q-bump code
      if residual^2>9*SigmaVtilde % bad voltage estimate by 3-SigmaX, bump Q 
        fprintf('Increasing SigmaX\n');
        SigmaX = SigmaX*2;
      end
      spkfData.SigmaX = SigmaX;      

      % Compute measurement-update estimate of SOC
      SOC = spkfData.SOC0 - spkfData.xhat(end)*(spkfData.Ts/(3600*spkfData.Q));         
  end
      
  % ----------------------------------------------------------------------
  % Steps 3a and 3b... operate only on four nearest neighbors
  Xind = getXind(Tk,SOC);
  X = getSigmaPoints(Xind);
  for theX = 1:size(X,2)
    [vhat,zk,Zsoc] = getVariables(X(:,theX),ik,Xind,Tk); 
    Z(:,theX) = [zk;vhat;Zsoc];
  end
  zk = Z*spkfData.alpham;
  Ztilde = zk(:,ones(1,size(Z,2))) - Z;
  SigmaZ = Ztilde * diag(spkfData.alphac) * Ztilde';
  boundzk = 3*sqrt(diag(SigmaZ));
    
  % ----------------------------------------------------------------------
  % Time to return 
  spkfData.priorI = ik;

  %% ======================================================================
  % The functions below this point implement the details of the higher-
  % level functionality indicated above
  % =======================================================================
  
  %% ----------------------------------------------------------------------
  % This function computes the indexing variables for the four blended mdls
  function Xind = getXind(Tk,SOC)
    % Find the two closest Zspts setpoints: "Zupper" and "Zlower"
    dZ = abs(SOC-Zpts); [~,iZ] = sort(dZ);
    Zupper = Zpts; iZupper = iZ; % default for single-model ROM
    Zlower = Zpts; iZlower = iZ; % default for single-model ROM
    if length(iZ)>1
      Zupper = Zpts(iZ(1)); iZupper = iZ(1);
      Zlower = Zpts(iZ(2)); iZlower = iZ(2);
      if Zupper<Zlower
        Zupper = Zpts(iZ(2)); iZupper = iZ(2);
        Zlower = Zpts(iZ(1)); iZlower = iZ(1);
      end
    end
    
    % Find the two closest temperature setpoints: "Tupper" and "Tlower"
    dT = abs(Tk-Tpts); [~,iT] = sort(dT);
    Tupper = Tpts; iTupper = iT; % default for single-model ROM
    Tlower = Tpts; iTlower = iT; % default for single-model ROM
    if length(iT)>1
      Tupper = Tpts(iT(1)); iTupper = iT(1);
      Tlower = Tpts(iT(2)); iTlower = iT(2);
      if Tupper<Tlower
        Tupper = Tpts(iT(2)); iTupper = iT(2);
        Tlower = Tpts(iT(1)); iTlower = iT(1);
      end
    end
    
    alphaZ = 0; alphaT = 0;
    if length(iZ)>1, alphaZ = (SOC-Zlower)/(Zupper-Zlower); end
    if length(iT)>1, alphaT = (Tk-Tlower)/(Tupper-Tlower); end    
    Xind.gamma = [(1-alphaT)*(1-alphaZ); % 00
                  (1-alphaT)*alphaZ;     % 01
                  alphaT*(1-alphaZ);     % 10
                  alphaT*alphaZ];        % 11
    Xind.theT = [iTlower iTlower iTupper iTupper]';
    Xind.theZ = [iZlower iZupper iZlower iZupper]';
  end

  %% ----------------------------------------------------------------------
  % This function computes sigma points, X, and delta between sigma points
  % and xhat, Xtilde for the four blended mdls
  function [X,Xtilde] = getSigmaPoints(Xind)
    indT = Xind.theT; indZ = Xind.theZ; 
    switch spkfData.method
      case 'OB'
        xh = [spkfData.M(indT(1),indZ(1)).xhat;
              spkfData.M(indT(2),indZ(2)).xhat;
              spkfData.M(indT(3),indZ(3)).xhat;
              spkfData.M(indT(4),indZ(4)).xhat;
              spkfData.x0];
        % faster to do 4 mchol small-matrix operations than one large one
        sigX = blkdiag(mchol(spkfData.M(indT(1),indZ(1)).SigmaX), ...
                       mchol(spkfData.M(indT(2),indZ(2)).SigmaX), ...
                       mchol(spkfData.M(indT(3),indZ(3)).SigmaX), ...
                       mchol(spkfData.M(indT(4),indZ(4)).SigmaX), ...
                       sqrt(spkfData.SigmaX0));
        X = xh(:,ones([1 8*n+3])) + spkfData.h*[zeros([4*n+1 1]), sigX, -sigX];
        Xtilde = xh(:,ones([1 8*n+3])) - X;
      case 'MB'
        xh = spkfData.xhat; sigX = mchol(spkfData.SigmaX);
        X = xh(:,ones([1 2*n+3])) + spkfData.h*[zeros([n+1 1]), sigX, -sigX];
        Xtilde = xh(:,ones([1 2*n+3])) - X;
    end
  end

  %% ----------------------------------------------------------------------
  % This function computes all cell nonlinear output variables; one set for
  % every input sigma point in X
% GLP NEEDS OB/MB variants
  function [Vcell,Z,Zsoc] = getVariables(X,ik,Xind,T)

    % ---------------------------------------------------------------------
    % Step 1: Extract the state sigma points of 4 nearest-neighbor models
    %         (needed for OB only)
    % ---------------------------------------------------------------------
    switch spkfData.method
      case 'OB'
        xk1 = X(0*n+1:1*n,:); % top n rows
        xk2 = X(1*n+1:2*n,:); % next n rows
        xk3 = X(2*n+1:3*n,:); % next n rows
        xk4 = X(3*n+1:4*n,:); % next n rows    
    end
    x0  = X(end,:);       % final row
    % Compute SOC for every sigma point
    xSOC = spkfData.SOC0 - X(end,:)*(spkfData.Ts/(3600*spkfData.Q)); 

    % ---------------------------------------------------------------------
    % Step 2: Initialize some variables/constants we will need
    % ---------------------------------------------------------------------
    SOCpAvg = cellData.function.pos.soc(xSOC,T);

    if any(SOCpAvg < 0)
      shortWarn('SOCpAvg < 0'); SOCpAvg(SOCpAvg<0) = 1e-6;
    end
    if any(SOCpAvg > 0.998)
      shortWarn('SOCpAvg > 1'); SOCpAvg(SOCpAvg>0.998) = 0.998;
    end
    
    % ---------------------------------------------------------------------
    % Step 3: Find the linear outputs as
    %         Z[k] = C*x[k] + D*iapp[k]
    % ---------------------------------------------------------------------
    Z = zeros(spkfData.nz,size(X,2));
    C1 = spkfData.M(Xind.theT(1),Xind.theZ(1)).C;
    C2 = spkfData.M(Xind.theT(2),Xind.theZ(2)).C;
    C3 = spkfData.M(Xind.theT(3),Xind.theZ(3)).C;
    C4 = spkfData.M(Xind.theT(4),Xind.theZ(4)).C;
    D1term = spkfData.M(Xind.theT(1),Xind.theZ(1)).D*ik;
    D2term = spkfData.M(Xind.theT(2),Xind.theZ(2)).D*ik;
    D3term = spkfData.M(Xind.theT(3),Xind.theZ(3)).D*ik;
    D4term = spkfData.M(Xind.theT(4),Xind.theZ(4)).D*ik;
    for k = 1:size(X,2)
      switch spkfData.method
        case 'OB'
          Z(:,k) = [C1*xk1(:,k) + D1term, C2*xk2(:,k) + D2term, ...
                    C3*xk3(:,k) + D3term, C4*xk4(:,k) + D4term]*Xind.gamma;
        case 'MB'
          Z(:,k) = [C1*X(1:n,k) + D1term, C2*X(1:n,k) + D2term, ...
                    C3*X(1:n,k) + D3term, C4*X(1:n,k) + D4term]*Xind.gamma;
      end
    end
    
    % --------------------------------------------------------
    % Step 4: Apply nonlinear corrections to the linear output
    % --------------------------------------------------------
    % No corrections needed for ifdl, if, or idl
    % But, we need If at current collectors later on
    If0 = Z(ind.If0,:);
    If3 = Z(ind.If3,:);

    % Solid surface stoichiometries (thetass)
    % (we are guaranteed that ind.negThetass and ind.posThetass not empty)
    % Note that we are implementing integrator-removed Thetass, so we need
    % to add in SOCpAvg instead of SOC0p
 
    Z(ind.posThetass,:) = Z(ind.posThetass,:) + SOCpAvg;
    if any(Z(ind.posThetass,:) < 0)
      shortWarn('posThetass < 0'); 
      ZposThetass = Z(ind.posThetass,:);
      ZposThetass(ZposThetass<0) = 1e-6;
      Z(ind.posThetass,:) = ZposThetass;
    end
    if any(Z(ind.posThetass,:) > 0.998)
      shortWarn('posThetass > 1'); 
      ZposThetass = Z(ind.posThetass,:);
      ZposThetass(ZposThetass>0.998) = 0.998;
      Z(ind.posThetass,:) = ZposThetass;
    end

    % Solid-electrolyte potential difference (phise)
    % The linear output from Z is integrator-removed version 
    % (we are guaranteed that ind.negPhise not empty)
    UocppAvg = cellData.function.pos.Uocp(SOCpAvg,T);
    Z(ind.negPhise,:) = Z(ind.negPhise,:);
    if ~isempty(ind.posPhise)
      Z(ind.posPhise,:) = Z(ind.posPhise,:) + UocppAvg;
    end
    % ### here
    % Compute electrolyte potential: first phie(0,t) then phie(1:3,t)
    % (we are guaranteed that ind.Phie not empty)
    PhieTilde3 = Z(ind.Phie(end),:);
    Phise0 = Z(ind.Phise0,:); 
    for theLoc = 1:length(ind.Phie)
      if loc.Phie(theLoc) ~= 0
        Z(ind.Phie(theLoc),:) = Z(ind.Phie(theLoc),:) - Phise0;
      end
    end

    % Compute electrolyte stoichiometries (thetae)
    % (we are guaranteed that ind.Thetae not empty)
    Z(ind.Thetae,:) = Z(ind.Thetae,:) + 1;
    if any(Z(ind.Thetae,:) < 0)
      shortWarn('Thetae < 0'); 
      ZThetae = Z(ind.Thetae,:);
      ZThetae(Zthetae < 0) = 1e-6;
      Z(ind,Thetae,:) = ZThetae;
    end

    % Compute overpotential at current-collectors via asinh method (eta)
    k0p = cellData.function.pos.k0(SOCpAvg,T);
    i0p = k0p*sqrt(Z(ind.Thetae(end),:).*...
               (1-Z(ind.Thetass3,:)).*Z(ind.Thetass3,:));
    
    Rfn    = cellData.function.neg.Rf(0.5,T);
    Rfp    = cellData.function.pos.Rf(SOCpAvg,T);
    negEta0 = -Z(ind.negPhie) - Rfn*ik;
    posEta3 = 2*R*T/F*asinh(If3./(2*i0p));

    % Compute cell voltage (ROMout.Vcell)
    Uocpp3 = cellData.function.pos.Uocp(Z(ind.Thetass3,:),T);
    

    Vcell = posEta3 - negEta0 + PhieTilde3 + Uocpp3 ...
            + (Rfp*Z(ind.Ifdl3,:) - Rfn*Z(ind.Ifdl0,:));

    % Compute solid potential (phis)
    Z(ind.posPhis,:) = Z(ind.posPhis,:) + Vcell;
    
    % Finally, compute SOC sigma points
    Zsoc = spkfData.SOC0 - x0*(spkfData.Ts/(3600*spkfData.Q));     
  end

  %% -----------------------------------------------------------------------
  % This function check the blending method the user has selected
  function checkBlendMethod
    % check to see what method is requested
    if strcmpi(spkfData.method,'MB')
%       error('Model-blended SPKF is not implemented yet.');
    end
  end

  %% ----------------------------------------------------------------------
  % This function sets up indices (ind) into the model linear output
  % vector "y" of variables needed to compute cell voltage and nonlinear
  % corrections. It also sets up locations (loc) in xtilde coordinates
  % of electrolyte variables.
  % -----------------------------------------------------------------------
  function setupIndsLocs
    % -- Find indices of outputs in model structure, to be used later 
    tfName = spkfData.tfData.names; % TF names
    tfLocs = spkfData.tfData.xLoc;  % TF regions (normalized x locations)

    % Negative electrode at x = 0
    ind.negIfdl    = find(strcmp(tfName,'negIfdl') == 1);
    ind.negIf      = find(strcmp(tfName,'negIf') == 1);
    ind.negPhise   = find(strcmp(tfName,'negPhise') == 1);
    loc.negIfdl    = tfLocs(ind.negIfdl);
    loc.negIf      = tfLocs(ind.negIf);
    loc.negPhise   = tfLocs(ind.negPhise);

    % Positive electrode
    ind.posIfdl    = find(strcmp(tfName,'posIfdl') == 1);
    ind.posIf      = find(strcmp(tfName,'posIf') == 1);
    ind.posPhis    = find(strcmp(tfName,'posPhis') == 1);
    ind.posPhise   = find(strcmp(tfName,'posPhise') == 1);
    ind.posThetass = find(strcmp(tfName,'posThetass') == 1);
    loc.posIfdl    = tfLocs(ind.posIfdl);
    loc.posIf      = tfLocs(ind.posIf);
    loc.posPhis    = tfLocs(ind.posPhis);
    loc.posPhise   = tfLocs(ind.posPhise);
    loc.posThetass = tfLocs(ind.posThetass);

    % Electrolyte potential across cell width
    ind.negPhie = find(strcmp(tfName,'negPhie') == 1);
    ind.dlPhie = find(strcmp(tfName,'dlPhie') == 1);
    ind.sepPhie = find(strcmp(tfName,'sepPhie') == 1); 
    ind.posPhie = find(strcmp(tfName,'posPhie') == 1);
    loc.negPhie = tfLocs(ind.negPhie);
    loc.dlPhie = tfLocs(ind.dlPhie);
    loc.sepPhie = tfLocs(ind.sepPhie); 
    loc.posPhie = tfLocs(ind.posPhie);

    % Electrolyte normalized concentration across cell width
    ind.dlThetae  = find(strcmp(tfName,'dlThetae')== 1);
    ind.sepThetae  = find(strcmp(tfName,'sepThetae')== 1); 
    ind.posThetae  = find(strcmp(tfName,'posThetae')== 1);
    loc.dlThetaes = tfLocs(ind.dlThetae);
    loc.sepThetaes = tfLocs(ind.sepThetae); 
    loc.posThetaes = tfLocs(ind.posThetae);

    % Combine inds and locs for variables across entire cell width
    ind.Ifdl    = [ind.negIfdl; ind.posIfdl];
    loc.Ifdl    = [loc.negIfdl; loc.posIfdl];
    ind.If      = [ind.negIf; ind.posIf];
    loc.If      = [loc.negIf; loc.posIf];
    ind.Phis    = [ind.posPhis];
    loc.Phis    = [loc.posPhis];
    ind.Phise   = [ind.negPhise; ind.posPhise];
    loc.Phise   = [loc.negPhise; loc.posPhise];
    ind.Thetass = [ind.posThetass];
    loc.Thetass = [loc.posThetass];
    
    ind.Phie    = [ind.negPhie;ind.dlPhie;ind.sepPhie;ind.posPhie];
    loc.Phie    = [loc.negPhie;loc.dlPhie;loc.sepPhie;loc.posPhie];
    ind.Thetae  = [ind.dlThetae;ind.sepThetae;ind.posThetae];
    loc.Thetae  = [loc.dlThetaes;loc.sepThetaes;loc.posThetaes];

    %-- Verify that variables required for nonlinear corrections exist 
    % Need to check:
    %  1. Ifdl    at both current-collectors
    %  2. If   at both current-collectors
    %  3. ROMout.Thetae  at both current-collectors
    %  4. Thetass at both current-collectors
    %  5. Phise   at negative-electrode current-collector 
    %  6. ROMout.Phie    at positive-electrode current-collector

    % Find location indexes
    ind.Ifdl0    = ind.Ifdl(loc.Ifdl == 0);
    ind.Ifdl3    = ind.Ifdl(loc.Ifdl == 3);
    ind.If0      = ind.If(loc.If == 0);
    ind.If3      = ind.If(loc.If == 3);
    ind.Thetass3 = ind.Thetass(loc.Thetass == 3);
    ind.Phise0   = ind.Phise(loc.Phise == 0);
    ind.Thetae0  = ind.Thetae(loc.Thetae == 0);
    ind.Thetae3  = ind.Thetae(abs(loc.Thetae - 3) < 0.001);

    % Check #1
    if isempty(ind.Ifdl0)
      error('Simulation requires ifdl at negative-collector!'); 
    end
    if isempty(ind.Ifdl3)
      error('Simulation requires ifdl at positive-collector!'); 
    end

    % Check #2
    if isempty(ind.If0)
      error('Simulation requires if at negative-collector!'); 
    end
    if isempty(ind.If3)
      error('Simulation requires if at positive-collector!'); 
    end

    % Check #3
    if loc.Thetae(1) > 0
      error('Simulation requires thetae at negative-collector!'); end
    if or(loc.Thetae(end)>3+eps,loc.Thetae(end)<3-eps)
      error('Simulation requires thetae at positive-collector!'); end

    % Check #4
    if isempty(ind.Thetass3)
      error('Simulation requires thetass at positive-collector!'); 
    end

    % Check #5
    if isempty(ind.Phise0)
      error('Simulation requires phise at negative-collector!'); 
    end

    % Check #6
    if loc.Phie(1) ~= 0
      error('Simulation requires phie at anode (x=0)!');

    end
    if or(loc.Phie(end)>3+eps,loc.Phie(end)<3-eps)
      error('Simulation requires phie at positive-collector!');
    end
  end
 
  %% ----------------------------------------------------------------------
  % This function displays a short warning message to the user if warnState
  % is []. (Does not display line numbers of warning, etc.)
  % -----------------------------------------------------------------------
  function shortWarn(msg)
    persistent warnState
    if strcmpi(msg,'on')
      warnState = []; 
    elseif strcmpi(msg,'off')
      warnState = 1;
    elseif isempty(warnState)
      cprintf([1,1/2,0],[' - Warning (' num2str(warnCount) '): ' msg '\n']);
      warnCount = warnCount + 1;
      if warnCount > 100
        warnState = 1;
      end
    end
  end

  %% -----------------------------------------------------------------------
  % M = mchol(G):  Given a symmetric matrix G, find a matrix E of "small" 
  % norm and L and D such that  G+E is positive definite, and G+E = L*D*L'.
  % Reference: Gill, Murray, and Wright, "Practical Optimization", p111.
  % Author: Brian Borchers (borchers@nmt.edu); speed enhancements by
  % Michael Zibulevsky. From "ldlt" package.
  function M = mchol(G)
    nn=size(G,1);        % size of matrix
    gamma=max(diag(G)); % compute quantities used by algorithm...
    zi=max(max(G-diag(diag(G))));
    nu=max([1,sqrt(nn^2-1)]);
    beta2=max([gamma, zi/nu, 1.0E-15]);

    CC=diag(diag(G)); % Init diag(CC) to diag(G)
    LL=zeros(nn); DD=zeros(nn); EE=zeros(nn);

    for j=1:nn % loop through, calculating column j of LL for d=1:nn
      bb=1:j-1; ee=j+1:nn;
      if (j > 1) % calculate jth row of LL
        LL(j,bb)=CC(j,bb)./diag(DD(bb,bb))';
      end
      if (j >= 2) % update jth column of CC
        if (j < nn) 
          CC(ee,j)=G(ee,j)-(LL(j,bb)*CC(ee,bb)')';
        end
      else
        CC(ee,j)=G(ee,j);
      end
      if (j == nn) % update theta
        theta=0;
      else
        theta=max(abs(CC(ee,j)));
      end
      DD(j,j)=max([eps;abs(CC(j,j));theta^2/beta2]); % update DD
      EE(j,j)=DD(j,j)-CC(j,j); % update EE
      indChol=(j*(nn+1)+1 : nn+1 : nn*nn)'; % update CC again
      CC(indChol)=CC(indChol)-(1/DD(j,j))*CC(ee,j).^2;
    end
    LL(1:nn+1:nn*nn)=1; % set diagonal of LL to 1s
    M = LL*sqrt(DD);
  end
end