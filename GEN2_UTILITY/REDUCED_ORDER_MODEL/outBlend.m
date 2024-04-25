% function ROMout = outBlend(simData,ROM,varargin)
% 
% Inputs:
%   simData  = input profile
%   ROM      = xRA pre-computed ROM model
%   Vmin     = optional: minimum permitted cell voltage
%   Vmax     = optional: maximum permitted cell voltage
%   Pmode    = optional: if set to 1, input is power and not current
%
% Outputs:
%   ROMout   = output structure contains simulated info
%
% This function simulates the physics-based ROM, using output-blending
%
% -- Changelog --
% 2023.09.30 | Compute i0p from MSMR | Wesley Hileman <whileman@uccs.edu>


function ROMout = outBlend(simData,ROM,varargin)
  Vmax = Inf; Vmin = -Inf; Pmode = 0;
  if nargin >= 4, Vmin = varargin{1}; Vmax = varargin{2}; end
  if nargin >=5
    Pmode = varargin{3}; 
    Vcell = (Vmin+Vmax)/2;
    options = optimoptions('fmincon','Display','off');
  end

  % Set some more accessible variables from the "ROM" data structure
  ROMmdls = ROM.ROMmdls; % precomputed state-space mdls at all setpoints
  cellData = ROM.cellData; % cell parameters

  % Re-define current and temperature data to enforce fixed Ts 
  [tk,ik,pk,Tk] = resampleInput;

  % Set up indices of simulation variables in linear outputs "y" 
  setupIndsLocs;
  
  % Create storage for simulation variables in output data structure
  setupOutputs;
  
  % Set up initial electrode local SOC (convert from percent)
  SOC0p = cellData.function.pos.soc(simData.SOC0/100,Tk(1));

  % Set up blending method and associated matrices/vectors
  [bigA,bigX,Tspts,Zspts,ZZ] = setupBlend;

  % Set up the cell state
  cellState.bigX = bigX;
  cellState.SOCpAvg = SOC0p;
  
  % Main simulation loop
  fprintf('Simulating ROM via outBlend...\n'); 
  for k = 0:length(ik)-1
    if Pmode % input to simulation is applied power
      shortWarn('off'); % no ROM warnings while seeking power level
      if pk(k+1) > 0 % search for ik to meet power spec...
        % initialize search using prior value of cell voltage
        ik(k+1) = fmincon(@(x) ...
                    (pk(k+1) - x*simStep(x,Tk(k+1),cellState))^2,...
                    pk(k+1)/Vcell,[],[],[],[],0,Inf,[],options);                          
      else        
        ik(k+1) = fmincon(@(x) ...
                    (pk(k+1) - x*simStep(x,Tk(k+1),cellState))^2,...
                    pk(k+1)/Vcell,[],[],[],[],-Inf,0,[],options);
      end
      shortWarn('on');
    end

    % first, attempt to update with user-supplied value of Iapp
    [Vcell,newCellState] = simStep(ik(k+1),Tk(k+1),cellState);

    % switch from current mode to constant-voltage mode
    if Vcell > Vmax
      shortWarn('off'); 
      ik(k+1) = fminbnd(@(x) (Vmax - simStep(x,Tk(k+1),cellState))^2,...
                  min(0,ik(k+1)),0);
      shortWarn('on'); 
      [Vcell,newCellState] = simStep(ik(k+1),Tk(k+1),cellState);
    end
    if Vcell < Vmin
      shortWarn('off'); 
      ik(k+1) = fminbnd(@(x) (Vmin - simStep(x,Tk(k+1),cellState))^2,...
                  0,max(0,ik(k+1))); 
      shortWarn('on'); 
      [Vcell,newCellState] = simStep(ik(k+1),Tk(k+1),cellState);
    end
    cellState = newCellState;

    % Update information to user every 100 iterations
    if (rem(k,100) == 0),fprintf('%d %2.8f\n',k,ROMout.cellSOC(k+1)); end    
  end
  
  % Return values as a structure
  storeData;

  % Finished! Now, return to the user
  
  %% ======================================================================
  % The functions below this point implement the details of the higher-
  % level functionality indicated above
  % =======================================================================

  %% ----------------------------------------------------------------------
  % This function resamples the simulation time/temperature/current vectors
  % to ensure an even sampling interval.
  % -----------------------------------------------------------------------
  function [tk,ik,pk,Tk] = resampleInput
    Ts = ROM.xraData.Tsamp; 
    tk = simData.time(1):Ts:simData.time(end); 
    if Pmode
      pk = interp1(simData.time,simData.Papp,tk); 
      ik = 0*pk;
    else
      ik = interp1(simData.time,simData.Iapp,tk); 
      pk = 0*ik;
    end
    Tk = interp1(simData.time,simData.T,tk)+273.15; 
    if size(ik,1) ~= size(Tk,1)
      error('Current and temperature profiles must be same length');
    end
  end

  %% ----------------------------------------------------------------------
  % This function sets up indices (ROM.ind) into the model linear output
  % vector "y" of variables needed to compute cell voltage and nonlinear
  % corrections. It also sets up locations (ROM.loc) in xtilde coordinates
  % of electrolyte variables.
  % -----------------------------------------------------------------------
  function setupIndsLocs
    % -- Find indices of outputs in model structure, to be used later 
    tfName = ROM.tfData.names; % TF names
    tfLocs = ROM.tfData.xLoc;  % TF regions (normalized x locations)
    ROM.tfLocs = tfLocs;

%     % Negative electrode at x = 0
    ROM.ind.negIfdl    = find(strcmp(tfName,'negIfdl') == 1);
    ROM.ind.negIf      = find(strcmp(tfName,'negIf') == 1);
%     ROM.ind.negIdl     = find(strcmp(tfName,'negIdl') == 1);
    ROM.ind.negPhise   = find(strcmp(tfName,'negPhise') == 1);


    % Positive electrode
    ROM.ind.posIfdl    = find(strcmp(tfName,'posIfdl') == 1);
    ROM.ind.posIf      = find(strcmp(tfName,'posIf') == 1);
    ROM.ind.posIdl     = find(strcmp(tfName,'posIdl') == 1);
    ROM.ind.posPhis    = find(strcmp(tfName,'posPhis') == 1);
    ROM.ind.posPhise   = find(strcmp(tfName,'posPhise') == 1);
    ROM.ind.posThetass = find(strcmp(tfName,'posThetass') == 1);

    % Electrolyte potential across cell width
    ROM.ind.negPhie = find(strcmp(tfName,'negPhie') == 1);
    ROM.ind.dllPhie = find(strcmp(tfName,'dllPhie') == 1);
    ROM.ind.sepPhie = find(strcmp(tfName,'sepPhie') == 1); 
    ROM.ind.posPhie = find(strcmp(tfName,'posPhie') == 1);
    ROM.loc.negPhie = tfLocs(ROM.ind.negPhie);
    ROM.loc.dllPhie = tfLocs(ROM.ind.dllPhie);
    ROM.loc.sepPhie = tfLocs(ROM.ind.sepPhie); 
    ROM.loc.posPhie = tfLocs(ROM.ind.posPhie);

    ROM.ind.Phie = [ROM.ind.negPhie; ROM.ind.dllPhie; ROM.ind.sepPhie; ROM.ind.posPhie];
    ROM.loc.Phie = [ROM.loc.negPhie; ROM.loc.dllPhie;ROM.loc.sepPhie;ROM.loc.posPhie];

    % Electrolyte normalized concentration across cell width
    ROM.ind.negThetae  = find(strcmp(tfName,'negThetae')== 1);
    ROM.ind.dllThetae  = find(strcmp(tfName,'dllThetae')== 1);
    ROM.ind.sepThetae  = find(strcmp(tfName,'sepThetae')== 1); 
    ROM.ind.posThetae  = find(strcmp(tfName,'posThetae')== 1);
    ROM.loc.negThetae = tfLocs(ROM.ind.negThetae);
    ROM.loc.dllThetae = tfLocs(ROM.ind.dllThetae);
    ROM.loc.sepThetae = tfLocs(ROM.ind.sepThetae); 
    ROM.loc.posThetae = tfLocs(ROM.ind.posThetae);

    ROM.ind.Thetae = [ROM.ind.negThetae;ROM.ind.dllThetae;ROM.ind.sepThetae;...
                       ROM.ind.posThetae];
    ROM.loc.Thetae = [ROM.loc.negThetae;ROM.loc.dllThetae;ROM.loc.sepThetae;...
                       ROM.loc.posThetae];

    %-- Verify that variables required for nonlinear corrections exist 
    % Need to check:
    %  1. Ifdl at both pos cc
    %  2. If at neg/pos cc
    %  3. Thetae at both cc
    %  4. Thetass at pos cc
    %  5. Phie at positive-electrode cc

    % Find location indexes
    ROM.ind.posIfdl3    = find(strcmp(tfName,'posIfdl')==1&tfLocs==3);
    ROM.ind.negIf0      = find(strcmp(tfName,'negIf')==1&tfLocs==0);
    ROM.ind.posIf3      = find(strcmp(tfName,'posIf')==1&tfLocs==3);
    ROM.ind.posThetass3 = find(strcmp(tfName,'posThetass')==1&tfLocs==3);
    ROM.ind.negThetae0  = find(strcmp(tfName,'negThetae')==1&tfLocs==0);
    ROM.ind.posThetae3  = find(strcmp(tfName,'posThetae')==1&tfLocs==3);
    ROM.ind.posPhie3    = find(strcmp(tfName,'posPhie')==1&tfLocs==3);

    % Check #1
    % NOTE: not required at neg for LMB!
    if isempty(ROM.ind.posIfdl3)
      error('Simulation requires ifdl at positive-collector!'); 
    end

    % Check #2
    if isempty(ROM.ind.negIf0)
      error('Simulation requires if at neg!'); 
    end
    if isempty(ROM.ind.posIf3)
      error('Simulation requires if at positive-collector!'); 
    end

    % Check #3
    if ROM.loc.Thetae(1) > 0
      error('Simulation requires thetae at neg!'); 
    end
    if or(ROM.loc.Thetae(end)>3+eps,ROM.loc.Thetae(end)<3-eps)
      error('Simulation requires thetae at positive-collector!'); 
    end

    % Check #4
    if isempty(ROM.ind.posThetass3)
      error('Simulation requires thetass at positive-collector!'); 
    end

    % Check #5
    if ROM.loc.Phie(1) == 0
      shortWarn('First phie x-location should not be zero. Ignoring');
      ROM.ind.Phie = ROM.ind.Phie(2:end);
    end
    if or(ROM.loc.Phie(end)>3+eps,ROM.loc.Phie(end)<3-eps)
      error('Simulation requires phie at positive-collector!');
    end
  end

  %% ----------------------------------------------------------------------
  % This function initializes storage for simulation variables in the
  % output data structure
  % -----------------------------------------------------------------------
  function setupOutputs
    duration = length(ik);
    
    % At negative electrode and its current collector
    ROMout.negIfdl    = zeros(duration,1);
    ROMout.negIf      = zeros(duration,length(ROM.ind.negIf));  
    ROMout.negPhie    = zeros(duration,length(ROM.ind.negPhie));
    ROMout.negPhise   = zeros(duration,length(ROM.ind.negPhise));

    % At positive electrode and its current collector
    ROMout.posIfdl    = zeros(duration,length(ROM.ind.posIfdl));
    ROMout.posIf      = zeros(duration,length(ROM.ind.posIf));
    ROMout.posIdl     = zeros(duration,length(ROM.ind.posIdl));
    ROMout.posPhis    = zeros(duration,length(ROM.ind.posPhis));
    ROMout.posPhise   = zeros(duration,length(ROM.ind.posPhise));
    ROMout.posThetass = zeros(duration,length(ROM.ind.posThetass));

    % Electrolyte variables
    ROMout.PhieTilde  = zeros(duration,length(ROM.ind.Phie));
    ROMout.Thetae     = zeros(duration,length(ROM.ind.Thetae));

    % Other cell and electrode quantities
    ROMout.Vcell   = zeros(duration,1);
    ROMout.cellSOC = zeros(duration,1);  
    ROMout.posSOC  = zeros(duration,1); 
  end      

  %% ----------------------------------------------------------------------
  % This function sets up the internal structure for performing the output
  % blending
  % -----------------------------------------------------------------------
  function [bigA,bigX,Tspts,Zspts,ZZ] = setupBlend
    % Extra all temperature and SOC setpoints (sort in ascending order)
    Tspts = sort(ROM.xraData.T)+273.15;
    Zspts = sort(ROM.xraData.SOC/100);

    % Create a "bigA" matrix.
    % Each column is the diagonal values of a specific A matrix.
    % The setpoint convention from column 1 to the last column is:
    % Z(1)T(1), Z(2)T(1), ..., Z(1)T(2), Z(2)T(2),..., Z(3)T(1),...
    % ROMmdls generated by xRA is a #Tspts by # Zspts structure
    [TT,ZZ] = size(ROMmdls);
    bigA = zeros(length(ROMmdls(1,1).A),TT*ZZ); 
    for tt = 1:TT
      for zz = 1:ZZ
        ind  = (tt-1)*ZZ + zz;
        ROMi = ROMmdls(tt,zz);
        bigA(:,ind) = real(diag(ROMi.A));

        % Don't forget to change res0 because we need only [phise]*
        ROMmdls(tt,zz).C(ROM.ind.posPhise,end) = 0;
      end
    end

    % Create the "bigX" matrix.
    % Each column stores the system states at a specific setpoint.
    bigX = zeros(size(bigA));
  end

  %% ----------------------------------------------------------------------
  % This function simulates one time step using output blending. Note that
  % while data are stored in the (k+1)st index, they can be overwritten if
  % the calling routine decides that an over/under current/power condition
  % has occured and re-calls this function
  % -----------------------------------------------------------------------
  function [Vcell,newCellState] = simStep(Iapp,T,cellState)
    F         = cellData.const.F;
    R         = cellData.const.R;    
    Q         = cellData.function.const.Q(); % cell capacity in Ah
    if isfield(cellData.function.const,'Rc')
        Rc = cellData.function.const.Rc();
    else
        Rc = 0;
    end
    theta0p   = cellData.function.pos.theta0();
    theta100p = cellData.function.pos.theta100();
    tauDLp    = cellData.function.pos.tauDL(SOC0p,T); 
    Cdlp      = cellData.function.pos.Cdl(SOC0p,T);
    nDLp      = cellData.function.pos.nDL();
    Cdldcp    = Cdlp^nDLp*tauDLp^(1-nDLp); % double-layer cap @ dc

    % Load present model state from "cellState"
    bigX      = cellState.bigX;
    SOCpAvg   = cellState.SOCpAvg;
    
    % ------------------------------------------
    % Step 1: Update cell SOC for next time step
    % ------------------------------------------
    % Store electrode average SOC data and compute cell SOC
    ROMout.posSOC(k+1)  = SOCpAvg;
    ROMout.cellSOC(k+1) = (SOCpAvg - theta0p)/(theta100p - theta0p);

    % Compute integrator input gains
    dUocppAvg = cellData.function.pos.dUocp(SOCpAvg,T);
    dQp       = abs(theta100p-theta0p);   
    res0p     = dQp/(3600*Q-Cdldcp*dQp*dUocppAvg);
   
    % Compute average positive-electrode SOC
    SOCpAvg = SOCpAvg + res0p*Iapp*(tk(2)-tk(1));
    if SOCpAvg < 0, shortWarn('Average SOCp < 0'); SOCpAvg = 0; end
    if SOCpAvg > 1, shortWarn('Average SOCp > 1'); SOCpAvg = 1; end

    % ---------------------------------------------------------------------
    % Step 2: Update linear output of closest four pre-computed ROM as
    %         y[k] = C*x[k] + D*iapp[k]
    % ---------------------------------------------------------------------
    % Find the two closest Zspts setpoints: "Zupper" and "Zlower"
    dZ = abs(ROMout.cellSOC(k+1)-Zspts); [~,iZ] = sort(dZ);
    Zupper = Zspts; iZupper = iZ;
    Zlower = Zspts; iZlower = iZ;
    if length(iZ)>1
      Zupper = Zspts(iZ(1)); iZupper = iZ(1);
      Zlower = Zspts(iZ(2)); iZlower = iZ(2);
      if Zupper<Zlower
        Zupper = Zspts(iZ(2)); iZupper = iZ(2);
        Zlower = Zspts(iZ(1)); iZlower = iZ(1);
      end
    end

    % Find the two closest temperature setpoints: "Tupper" and "Tlower"
    dT = abs(T-Tspts); [~,iT] = sort(dT);
    Tupper = Tspts; iTupper = iT;
    Tlower = Tspts; iTlower = iT;
    if length(iT)>1
      Tupper = Tspts(iT(1)); iTupper = iT(1);
      Tlower = Tspts(iT(2)); iTlower = iT(2);
      if Tupper<Tlower
        Tupper = Tspts(iT(2)); iTupper = iT(2);
        Tlower = Tspts(iT(1)); iTlower = iT(1);
      end
    end

    % Compute the four outputs at the time sample
    % Each yk is a column matrix (#TFlocs by 1)
    yk1 = ROMmdls(iTlower,iZlower).C*bigX(:,(iTlower-1)*ZZ+iZlower)...
          + ROMmdls(iTlower,iZlower).D*Iapp;
    yk2 = ROMmdls(iTlower,iZupper).C*bigX(:,(iTlower-1)*ZZ+iZupper)...
          + ROMmdls(iTlower,iZupper).D*Iapp;
    yk3 = ROMmdls(iTupper,iZlower).C*bigX(:,(iTupper-1)*ZZ+iZlower)...
          + ROMmdls(iTupper,iZlower).D*Iapp;
    yk4 = ROMmdls(iTupper,iZupper).C*bigX(:,(iTupper-1)*ZZ+iZupper)...
          + ROMmdls(iTupper,iZupper).D*Iapp;

    % ---------------------------------------------------
    % Step 3: Update states of all pre-computed ROM as
    %         bigX[k+1] = bigA.*bigX[k] + B*iapp[k]
    % ---------------------------------------------------
    bigX = bigA.*bigX + Iapp;
        
    % ----------------------------------------------------------
    % Step 4: Compute the linear output at the present setpoint
    % ----------------------------------------------------------
    % Compute the blending factors
    alphaZ = 0; alphaT = 0;
    if length(iZ)>1
      alphaZ = (ROMout.cellSOC(k+1)-Zlower)/(Zupper-Zlower); 
    end
    if length(iT)>1
      alphaT = (T-Tlower)/(Tupper-Tlower); 
    end

    % Compute the present linear output yk (a column vector)
    yk = (1-alphaT)*((1-alphaZ)*yk1 + alphaZ*yk2)...
         +alphaT*((1-alphaZ)*yk3 + alphaZ*yk4);

    % -------------------------------------------------------
    % Step 5: Add nonlinear corrections to the linear output
    % -------------------------------------------------------

    % Compute variables at specific x-locations.
    Thetae0 = yk(ROM.ind.negThetae0) + 1;
    Thetae3 = yk(ROM.ind.posThetae3) + 1;
    Thetass3 = yk(ROM.ind.posThetass3) + SOC0p;
    If0 = yk(ROM.ind.negIf0);
    If3 = yk(ROM.ind.posIf3);
    PhieTilde3 = yk(ROM.ind.posPhie3);
    Ifdl3 = yk(ROM.ind.posIfdl3);
    % Compute overpotential at current-collectors via asinh method (eta)
    if isfield(cellData.function.pos,'U0')
        % MSMR kinetics.
        electrode = MSMR(cellData.function.pos);
        ctData = electrode.Rct( ...
            cellData.function.pos,'theta',Thetass3);
        i0p = ctData.i0*sqrt(Thetae3);
    else
        k0p = cellData.function.pos.k0(ROMout.posSOC(k+1),T);
        i0p = k0p*sqrt(Thetae3*(1-Thetass3)*Thetass3);
    end
    k0n = cellData.function.neg.k0(ROMout.posSOC(k+1),T);
    i0n = k0n*sqrt(Thetae0);
    posEta3 = 2*R*T/F*asinh(If3/(2*i0p));
    negEta0 = 2*R*T/F*asinh(If0/(2*i0n));

    % Interfacial total molar rate (ifdl)
    ROMout.posIfdl(k+1,:) = yk(ROM.ind.posIfdl);

    % Interfacial faradaic molar rate (if)
    ROMout.negIf(k+1,:)   = yk(ROM.ind.negIf);
    ROMout.posIf(k+1,:)   = yk(ROM.ind.posIf);

    % Interfacial nonfaradaic molar rate (Idl)
%     ROMout.negIdl(k+1,:)  = yk(ROM.ind.negIdl);
%     ROMout.posIdl(k+1,:)  = yk(ROM.ind.posIdl);

    % Solid surface stoichiometries (thetass)
    ROMout.posThetass(k+1,:) = yk(ROM.ind.posThetass) + SOC0p;
    if any(ROMout.posThetass(k+1,:) < 0)
      shortWarn('posThetass < 0'); 
      ROMout.posThetass(ROMout.posThetass < 0) = 1e-6;
    end
    if any(ROMout.posThetass(k+1,:) > 1)
      shortWarn('posThetass > 1'); 
      ROMout.posThetass(ROMout.posThetass > 1) = 1-1e-6;
    end

    % Solid-electrolyte potential difference (phise)
    % The linear output from yk is integrator-removed version
    UocppAvg = cellData.function.pos.Uocp(ROMout.posSOC(k+1),T);
    ROMout.negPhise(k+1,:) = yk(ROM.ind.negPhise);
    ROMout.posPhise(k+1,:) = yk(ROM.ind.posPhise)  + UocppAvg;

    % Compute electrolyte potential.
    ROMout.PhieTilde(k+1,:) = yk(ROM.ind.Phie);

    % Compute electrolyte stoichiometries (thetae)
    ROMout.Thetae(k+1,:) = yk(ROM.ind.Thetae) + 1;
    if any(ROMout.Thetae(k+1,:) < 0)
      shortWarn('Thetae < 0'); ROMout.Thetae(ROMout.Thetae < 0) = 1e-6;
    end

    % Compute cell voltage (ROMout.Vcell)
    Uocpp3 = cellData.function.pos.Uocp(Thetass3,T);
    Rfn    = cellData.function.neg.Rf(0.5,T);
    Rfp    = cellData.function.pos.Rf(ROMout.posSOC(k+1),T);
    ROMout.Vcell(k+1) = ...
        (posEta3+Uocpp3+Rfp*Ifdl3)-(negEta0+Rfn*Iapp)+PhieTilde3;

    % Finally, compute solid potential (phis)
    ROMout.posPhis(k+1,:) = yk(ROM.ind.posPhis) + ROMout.Vcell(k+1);

    % Correct Vcell for tab resistance (must do after Phis computation!).
    ROMout.Vcell(k+1) =  ROMout.Vcell(k+1) - Iapp*Rc;
    
    % Function outputs: voltage and updated cell state
    Vcell = ROMout.Vcell(k+1);
    newCellState.bigX = bigX;
    newCellState.SOCpAvg = SOCpAvg;
  end

  %% ----------------------------------------------------------------------
  % This function packages all simulation data into the output data
  % structure to return to the user
  % -----------------------------------------------------------------------
  function storeData
    % A final information update
    fprintf('%d %2.8f \n Voltage: %2.8f\n',length(ik),...
            ROMout.cellSOC(end),ROMout.Vcell(end));

    % Save basic information
    ROMout.blending = 'output-blending';
    ROMout.cellData = cellData;
    ROMout.time     = tk(:);
    ROMout.Iapp     = ik(:);
    ROMout.T        = Tk(:);

    % Save each variable
    ROMout.Ifdl    = [ROMout.posIfdl];
    ROMout.If      = [ROMout.negIf ROMout.posIf];
    ROMout.Idl     = [ROMout.posIdl];
    ROMout.Phis    = [ROMout.posPhis];
    ROMout.Phise   = [ROMout.negPhise ROMout.posPhise];
    ROMout.Thetass = [ROMout.posThetass];

    % Save locations of each variable
    tfLocs               = ROM.tfLocs;
    ROMout.xLocs.Ifdl    = [tfLocs(ROM.ind.posIfdl)];
    ROMout.xLocs.If      = [tfLocs(ROM.ind.negIf); tfLocs(ROM.ind.posIf)];
    ROMout.xLocs.Phis    = [tfLocs(ROM.ind.posPhis)];
    ROMout.xLocs.PhieTilde = [tfLocs(ROM.ind.negPhie);
                            tfLocs(ROM.ind.dllPhie);...
                            tfLocs(ROM.ind.sepPhie);...
                            tfLocs(ROM.ind.posPhie)];
    ROMout.xLocs.Phise   = [tfLocs(ROM.ind.negPhise);...
                            tfLocs(ROM.ind.posPhise)];
    ROMout.xLocs.Phise   = [tfLocs(ROM.ind.negPhise) ; tfLocs(ROM.ind.posPhise)];
    ROMout.xLocs.Thetass = [tfLocs(ROM.ind.posThetass)];
    ROMout.xLocs.Thetae  = [tfLocs(ROM.ind.negThetae);
                            tfLocs(ROM.ind.dllThetae);...
                            tfLocs(ROM.ind.sepThetae);...
                            tfLocs(ROM.ind.posThetae)];
  end

  %% ----------------------------------------------------------------------
  % This function displays a short warning message to the user if warnState
  % is []. (Does not display line numbers of warning, etc.)
  % -----------------------------------------------------------------------
  function shortWarn(msg)
    if isfield(simData,'warnOff'), return; end
    persistent warnState
    if strcmpi(msg,'on')
      warnState = []; 
    elseif strcmpi(msg,'off')
      warnState = 1;
    elseif isempty(warnState)
      cprintf([1,1/2,0],[' - Warning: ' msg '\n']);
    end
  end
end