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
    ROM.ind.dlPhie = find(strcmp(tfName,'dlPhie') == 1);
    ROM.ind.sepPhie = find(strcmp(tfName,'sepPhie') == 1); 
    ROM.ind.posPhie = find(strcmp(tfName,'posPhie') == 1);
    ROM.loc.negPhie = tfLocs(ROM.ind.negPhie);
    ROM.loc.dlPhie = tfLocs(ROM.ind.dlPhie);
    ROM.loc.sepPhie = tfLocs(ROM.ind.sepPhie); 
    ROM.loc.posPhie = tfLocs(ROM.ind.posPhie);

    ROM.ind.Phie = [ROM.ind.negPhie; ROM.ind.dlPhie;ROM.ind.sepPhie;ROM.ind.posPhie];
    ROM.loc.Phie = [ROM.loc.negPhie; ROM.loc.dlPhie;ROM.loc.sepPhie;ROM.loc.posPhie];

    % Electrolyte normalized concentration across cell width
    ROM.ind.dlThetae  = find(strcmp(tfName,'dlThetae')== 1);
    ROM.ind.sepThetae  = find(strcmp(tfName,'sepThetae')== 1); 
    ROM.ind.posThetae  = find(strcmp(tfName,'posThetae')== 1);
    ROM.loc.dlThetae = tfLocs(ROM.ind.dlThetae);
    ROM.loc.sepThetae = tfLocs(ROM.ind.sepThetae); 
    ROM.loc.posThetae = tfLocs(ROM.ind.posThetae);

    ROM.ind.Thetae = [ROM.ind.dlThetae;ROM.ind.sepThetae;...
                       ROM.ind.posThetae];
    ROM.loc.Thetae = [ROM.loc.dlThetae;ROM.loc.sepThetae;...
                       ROM.loc.posThetae];

    %-- Verify that variables required for nonlinear corrections exist 
    % Need to check:
    %  1. Ifdl    at both current-collectors
    %  2. If   at both current-collectors
    %  3. ROMout.Thetae  at both current-collectors
    %  4. Thetass at both current-collectors
    %  5. Phise   at negative-electrode current-collector 
    %  6. ROMout.Phie    at positive-electrode current-collector

    % Find location indexes
    ROM.ind.negIfdl0    = find(strcmp(tfName,'negIfdl') == 1 ...
                                & tfLocs == 0);
    ROM.ind.posIfdl3    = find(strcmp(tfName,'posIfdl') == 1 ...
                                & tfLocs == 3);
    ROM.ind.negIf0   = find(strcmp(tfName,'negIf') == 1);
    ROM.ind.posIf3   = find(strcmp(tfName,'posIf') == 1 ...
                                & tfLocs == 3);
    ROM.ind.posThetass3 = find(strcmp(tfName,'posThetass') == 1 ...
                                & tfLocs == 3);
    ROM.ind.negPhise0   = find(strcmp(tfName,'negPhise') == 1 ...
                                & tfLocs == 0);

    % Check #1
%     if isempty(ROM.ind.negIfdl0)
%       error('Simulation requires ifdl at negative-collector!'); 
%     end
    if isempty(ROM.ind.posIfdl3)
      error('Simulation requires ifdl at positive-collector!'); 
    end

    % Check #2
%     if isempty(ROM.ind.negIf0)
%       error('Simulation requires if at negative-collector!'); 
%     end
    if isempty(ROM.ind.posIf3)
      error('Simulation requires if at positive-collector!'); 
    end

    % Check #3
    if ROM.loc.Thetae(1) > 0
      error('Simulation requires thetae at negative-collector!'); end
    if or(ROM.loc.Thetae(end)>3+eps,ROM.loc.Thetae(end)<3-eps)
      error('Simulation requires thetae at positive-collector!'); end

    % Check #4
    if isempty(ROM.ind.posThetass3)
      error('Simulation requires thetass at positive-collector!'); 
    end

    % Check #5
%     if isempty(ROM.ind.negPhise0)
%       error('Simulation requires phise at negative-collector!'); 
%     end

    % Check #6
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
    ROMout.negIfdl    = zeros(duration,length(ROM.ind.negIfdl));
    ROMout.negIf      = zeros(duration,length(ROM.ind.negIf));
%     ROMout.negIdl     = zeros(duration,length(ROM.ind.negIdl));   
    ROMout.negPhie   = zeros(duration,length(ROM.ind.negPhie));
    ROMout.negPhise   = zeros(duration,length(ROM.ind.negPhise));
    ROMout.negIfdl0    = zeros(duration,length(ROM.ind.negIfdl0));
    ROMout.negIf0   = zeros(duration,length(ROM.ind.negIf0));
%     ROMout.negEta0     = zeros(duration,length(ROM.ind.negIf0));
    ROMout.negPhise0   = zeros(duration,length(ROM.ind.negPhise0));


    % At positive electrode and its current collector
    ROMout.posIfdl    = zeros(duration,length(ROM.ind.posIfdl));
    ROMout.posIf      = zeros(duration,length(ROM.ind.posIf));
    ROMout.posIdl     = zeros(duration,length(ROM.ind.posIdl));
    ROMout.posPhis    = zeros(duration,length(ROM.ind.posPhis));
    ROMout.posPhise   = zeros(duration,length(ROM.ind.posPhise));
    ROMout.posThetass = zeros(duration,length(ROM.ind.posThetass));

    ROMout.posIfdl3    = zeros(duration,length(ROM.ind.posIfdl3));
    ROMout.posIf3   = zeros(duration,length(ROM.ind.posIf3));
    ROMout.posEta3     = zeros(duration,length(ROM.ind.posIf3));
    ROMout.posThetass3 = zeros(duration,length(ROM.ind.posThetass3));

    % Electrolyte variables
    ROMout.Phie   = zeros(duration,length(ROM.ind.Phie)+1); % add x=0 loc
    ROMout.Thetae = zeros(duration,length(ROM.ind.Thetae));

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
%         ROMmdls(tt,zz).C(ROM.ind.negPhise,end) = 0;
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
        Rc    = cellData.function.const.Rc();
    else
        Rc    = 0;
    end
    theta0p   = cellData.function.pos.theta0();
    theta100p = cellData.function.pos.theta100();
    
%     wDLn      = cellData.function.neg.wDL(SOC0n,T); 
    wDLp      = cellData.function.pos.wDL(SOC0p,T); 
%     Cdln      = cellData.function.neg.Cdl(SOC0n,T); 
    Cdlp      = cellData.function.pos.Cdl(SOC0p,T);
%     nDLn      = cellData.function.neg.nDL();
    nDLp      = cellData.function.pos.nDL();
%     Cdleffn   = (Cdln^(2-nDLn))*(wDLn^(nDLn-1));
    Cdleffp   = (Cdlp^(2-nDLp))*(wDLp^(nDLp-1));

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
    res0p     =  dQp/(3600*Q-Cdleffp*dQp*dUocppAvg);
   
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
    % Interfacial total molar rate (ifdl)
    ROMout.negIfdl(k+1,:) = yk(ROM.ind.negIfdl);
    ROMout.posIfdl(k+1,:) = yk(ROM.ind.posIfdl);
    ROMout.negIfdl0(k+1)  = yk(ROM.ind.negIfdl0);
    ROMout.posIfdl3(k+1)  = yk(ROM.ind.posIfdl3);

    % Interfacial faradaic molar rate (if)
%     ROMout.negIf(k+1,:)   = yk(ROM.ind.negIf);
    ROMout.posIf(k+1,:)   = yk(ROM.ind.posIf);
    ROMout.negIf0(k+1)    = yk(ROM.ind.negIf0);
    ROMout.posIf3(k+1)    = yk(ROM.ind.posIf3);

    % Interfacial nonfaradaic molar rate (Idl)
%     ROMout.negIdl(k+1,:)  = yk(ROM.ind.negIdl);
    ROMout.posIdl(k+1,:)  = yk(ROM.ind.posIdl);

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


    ROMout.posThetass3(k+1) = yk(ROM.ind.posThetass3) + SOC0p;
    if any(ROMout.posThetass3(k+1) < 0)
      shortWarn('posThetass3 < 0'); 
      ROMout.posThetass3(ROMout.posThetass3 < 0) = 1e-6;
    end
    if any(ROMout.posThetass3(k+1) > 1)
      shortWarn('posThetass3 > 1'); 
      ROMout.posThetass3(ROMout.posThetass3 > 1) = 1-1e-6;
    end

    % Solid-electrolyte potential difference (phise)
    % The linear output from yk is integrator-removed version
    UocppAvg = cellData.function.pos.Uocp(ROMout.posSOC(k+1),T);
%     ROMout.negPhise(k+1,:) = yk(ROM.ind.negPhise)  + UocpnAvg;
    ROMout.posPhise(k+1,:) = yk(ROM.ind.posPhise)  + UocppAvg;
    ROMout.negPhise0(k+1)  = yk(ROM.ind.negPhise0) ;

    % Compute electrolyte potential: first phie(0,t) then phie(1:3,t)
    ROMout.Phie(k+1,1)     = yk(ROM.ind.negPhie); 
    ROMout.Phie(k+1,2:end) = yk(ROM.ind.Phie) + yk(ROM.ind.negPhie); 

    % Compute electrolyte stoichiometries (thetae)
    ROMout.Thetae(k+1,:) = yk(ROM.ind.Thetae) + 1;
    if any(ROMout.Thetae(k+1,:) < 0)
      shortWarn('Thetae < 0'); ROMout.Thetae(ROMout.Thetae < 0) = 1e-6;
    end

    % Compute overpotential at current-collectors via asinh method (eta)
    k0p = cellData.function.pos.k0(ROMout.posSOC(k+1),T);
    if length(k0p)>1
        % MSMR kinetics.
        electrode = MSMR(cellData.function.pos);
        ctData = electrode.Rct(cellData.function.pos,'theta',SOCpAvg);
        i0p = ctData.i0;
    else
        i0p = k0p*sqrt(ROMout.Thetae(k+1,end)*...
                   (1-ROMout.posThetass3(k+1))*ROMout.posThetass3(k+1));
    end
    
    ROMout.posEta3(k+1) = 2*R*T/F*asinh(ROMout.posIf3(k+1)/(2*i0p));

    % Compute cell voltage (ROMout.Vcell)
    Uocpp3 = cellData.function.pos.Uocp(ROMout.posThetass3(k+1),T);
    Rfn    = cellData.function.neg.Rf(0.5,T);
    ROMout.negEta0(k+1) = -ROMout.Phie(k+1,1) - Rfn*Iapp;
    Rfp    = cellData.function.pos.Rf(ROMout.posSOC(k+1),T);

    ROMout.Vcell(k+1) = ROMout.posEta3(k+1) - ROMout.negEta0(k+1)...
        + yk(ROM.ind.Phie(end)) + Uocpp3 ...
        + (Rfp*ROMout.posIfdl3(k+1) - Rfn*ROMout.negIfdl(k+1));

    % Finally, compute solid potential (phis)
    ROMout.posPhis(k+1,:) = yk(ROM.ind.posPhis) + ROMout.Vcell(k+1);
    
    % Update Vcell if Rc is nonzero
    ROMout.Vcell(k+1) = ROMout.Vcell(k+1) - Rc*Iapp;
    
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
    ROMout.Ifdl    = [ROMout.negIfdl0 ROMout.posIfdl];
%     ROMout.Ifdl    = [ROMout.negIfdl ROMout.posIfdl];
    ROMout.If      = [ ROMout.negIf0 ROMout.posIf]; % [ROMout.negIf ROMout.posIf];
    ROMout.Idl     = [ ROMout.posIdl];% [ROMout.negIdl ROMout.posIdl];
    ROMout.Phis    = [ROMout.posPhis];
    ROMout.Phise   = [ROMout.negPhise0 ROMout.posPhise]; % [ROMout.negPhise ROMout.posPhise];
    ROMout.Thetass = [ROMout.posThetass];

    % Save locations of each variable
    tfLocs               = ROM.tfLocs;
%     ROMout.xLocs.Ifdl    = [tfLocs(ROM.ind.negIfdl);...
%                             tfLocs(ROM.ind.posIfdl)];
%     ROMout.xLocs.If      = [tfLocs(ROM.ind.negIf);...
%                             tfLocs(ROM.ind.posIf)];
%     ROMout.xLocs.Idl     = [tfLocs(ROM.ind.negIdl);...
%                             tfLocs(ROM.ind.posIdl)];
    ROMout.xLocs.Ifdl    = [tfLocs(ROM.ind.negIfdl); tfLocs(ROM.ind.posIfdl)];
    ROMout.xLocs.If      = [tfLocs(ROM.ind.negIf); tfLocs(ROM.ind.posIf)];
%     ROMout.xLocs.Idl     = [tfLocs(ROM.ind.posIdl)];
    ROMout.xLocs.Phis    = [tfLocs(ROM.ind.posPhis)];
    ROMout.xLocs.Phie    = [tfLocs(ROM.ind.negPhie);tfLocs(ROM.ind.dlPhie);...
                            tfLocs(ROM.ind.sepPhie);...
                            tfLocs(ROM.ind.posPhie)];
    ROMout.xLocs.Phise   = [tfLocs(ROM.ind.negPhise);...
                            tfLocs(ROM.ind.posPhise)];
   ROMout.xLocs.Phise   = [tfLocs(ROM.ind.negPhise) ; tfLocs(ROM.ind.posPhise)];
    ROMout.xLocs.Thetass = [tfLocs(ROM.ind.posThetass)];
    ROMout.xLocs.Thetae  = [tfLocs(ROM.ind.dlThetae);...
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