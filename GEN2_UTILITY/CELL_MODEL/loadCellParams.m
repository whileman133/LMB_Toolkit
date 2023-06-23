% function cellParams = loadCellParams(fileName)
% 
% Input:
%   fileName   = full path to Excel spreadsheet of cell parameter values
% Output:
%   cellParams = data structure containing all parameter values

%
% This file is provided as a supplement to: Plett, Gregory L. and Trimboli,
% M. Scott, "Battery Management Systems, Volume III, Physics-Based
% Methods," Artech House, 2021. 
%
% -- Changelog --
% 07.06.2022 | Wesley Hileman | Edited, Modified for LMB
%   - readParamTable: use `readcell` instead of `xlsread` (fixes sheet 
%     detection issue when reading model created by saveCellParams)
%   - loadTable: use `readmatrix` instead of `xlsread` (fixes sheet 
%     detection issue when reading model created by saveCellParams)
%   - Extract part of readParamTable into externally-accessible 
%     loadParamTable function.
%   - convertStandardToLumped: 
%       * Remove Rct computation for negative electrode; already computed 
%         in tfCommon. Remove nF from negative electrode (no solid-phase 
%         diffusion at the anode).
%       * Remove temperature dependence from psi (always appears as psi*T
%         in the model).
% -- Changelog --
% 07.19.2022 | Aloisio | 
%   - changed the SOCvec to create LUT of OCP 0.0001

function cellParams = loadCellParams(fileName)
  inputCellData = readParamTable(fileName,'Parameters'); % read Excel sheet
  if ~inputCellData.lumped % standard parameters, convert to lumped
    inputCellData = convertStandardToLumped(inputCellData);
  end
  cellParams = cellStr2Fun(inputCellData,1); % convert strings to functions
  if cellParams.MSMR
    cellParams.function.pos.Uocp  = makeOCP(cellParams.function.pos);
    cellParams.function.pos.dUocp = makedOCP(cellParams.function.pos);
  end
  cellParams = makeOCV(cellParams);
end

% Read the Excel spreadsheet from file=fileName, tab=sheet
function cellParams = readParamTable(fileName,sheet)
  cellParams = []; % initialize to empty structure
  
  % Read contents of Worksheet into cell array.
  [ind,data] = loadParamTable(fileName,sheet);
  
  % Individually process the main sections of the XLSX data file
  cellParams.const.F = 96485.3365; % Faraday's constant
  cellParams.const.R = 8.3144621;  % Universal gas constant
  cellParams.const.T = 298.15;
  cellParams = processGen(data(ind.gen,:),cellParams);
  cellParams = processCst(data(ind.const,:),cellParams,fileName);
  cellParams = processNeg(data(ind.neg,:),cellParams,fileName);
  cellParams = processDL(data(ind.DL,:),cellParams,fileName);
  cellParams = processSep(data(ind.sep,:),cellParams,fileName);
  cellParams = processPos(data(ind.pos,:),cellParams,fileName);

  % Add SOC function to each electrode
  theta0_pos = str2func(cellParams.function.pos.theta0); 
  theta100_pos = str2func(cellParams.function.pos.theta100);
  cellParams.function.pos.soc = sprintf('@(x,T) (%g + x*(%g))',...
    theta0_pos(),theta100_pos()-theta0_pos());  
end

% Process the "#general" segment of the data file
% ... cell name and type (lumped or standard)
function cellParams = processGen(data,cellParams)
  % Prepare data by extracting only relevant part, deleting "Code Name"
  vars = data(:,2); vals = data(:,3);
  indVar = cellfun(@isstr,vars); % index of strings in "vars" column
  indVar2 = logical(indVar.*indVar([end,1:end-1])); % Delete "Code Name"
  vars = vars(indVar2); vals = vals(indVar2); 
  indName = strcmpi(vars,'name'); % Search for cell-name row
  if sum(indName) == 0
    error('Missing "name" entry in the general section of xlsx file');
  end
  cellParams.name = vals{indName};
  indLumped = strcmpi(vars,'lumped'); % Search for cell-type row
  if isempty(indLumped)
    error('Missing "lumped" entry in the general section of xlsx file');
  end
  cellParams.lumped = vals{indLumped};
  indMSMR = strcmpi(vars,'MSMR'); % Search for kinetics-type row
  if sum(indMSMR) == 0
    cellParams.MSMR = 0;
  else
    cellParams.MSMR = vals{indMSMR};
  end
  indTM = strcmpi(vars,'TM'); % Search for kinetics-type row
  if sum(indTM) == 0
    cellParams.TM = 0;
  else
    cellParams.TM = vals{indTM};
  end
end

% Process the "#const" segment of the data file
function cellParams = processCst(data,cellParams,fileName)
  [vars,vals,Eact] = parseData(data,fileName);
  for j = 1:length(vars)
    const.(vars{j}) = vals{j};
    if Eact{j} ~= 0 && ~isnan(Eact{j})
      const.(sprintf('%s__Ea',vars{j})) = sprintf('%g',1000*Eact{j});
    end
  end
  cellParams.function.const = const;
%   cellParams.function.const = addTemps(const);
end

% Process the "#neg" segment of the data file
function cellParams = processNeg(data,cellParams,fileName)
  [vars,vals,Eact] = parseData(data,fileName);
  for j = 1:length(vars)
    neg.(vars{j}) = vals{j};
    if ischar(Eact{j}) % then probably a vector
      EactVect = str2num(Eact{j}); %#ok<ST2NM>
      EactVect = 1000*EactVect(:);
      EactStr = sprintf(' %g;',EactVect);
      neg.(sprintf('%s__Ea',vars{j})) = sprintf('[%s]',EactStr(1:end-1)); % delete ";" from end
    elseif Eact{j} ~= 0 && ~isnan(Eact{j})
      neg.(sprintf('%s__Ea',vars{j})) = sprintf('%g',1000*Eact{j});
    end
  end
  cellParams.function.neg = neg;
end

% Process the "#DL" segment of the data file
function cellParams = processDL(data,cellParams,fileName)
  [vars,vals,Eact] = parseData(data,fileName);
  for j = 1:length(vars)
    DL.(vars{j}) = vals{j};
    if ischar(Eact{j}) % then probably a vector
      EactVect = str2num(Eact{j}); %#ok<ST2NM>
      EactVect = 1000*EactVect(:);
      EactStr = sprintf(' %g;',EactVect);
      DL.(sprintf('%s__Ea',vars{j})) = sprintf('[%s]',EactStr(1:end-1)); % delete ";" from end
    elseif Eact{j} ~= 0 && ~isnan(Eact{j})
      DL.(sprintf('%s__Ea',vars{j})) = sprintf('%g',1000*Eact{j});
    end
  end
  cellParams.function.DL = DL;
end

% Process the "#sep" segment of the data file
function cellParams = processSep(data,cellParams,fileName)
  [vars,vals,Eact] = parseData(data,fileName);
  for j = 1:length(vars)
    sep.(vars{j}) = vals{j};
    if ischar(Eact{j}) % then probably a vector
      EactVect = str2num(Eact{j}); %#ok<ST2NM>
      EactVect = 1000*EactVect(:);
      EactStr = sprintf(' %g;',EactVect);
      sep.(sprintf('%s__Ea',vars{j})) = sprintf('[%s]',EactStr(1:end-1)); % delete ";" from end
    elseif Eact{j} ~= 0 && ~isnan(Eact{j})
      sep.(sprintf('%s__Ea',vars{j})) = sprintf('%g',1000*Eact{j});
    end
  end
  cellParams.function.sep = sep;
end

% Process the "#pos" segment of the data file
function cellParams = processPos(data,cellParams,fileName)
  [vars,vals,Eact] = parseData(data,fileName);
  for j = 1:length(vars)
    pos.(vars{j}) = vals{j};
    if ischar(Eact{j}) % then probably a vector
      EactVect = str2num(Eact{j}); %#ok<ST2NM>
      EactVect = 1000*EactVect(:);
      EactStr = sprintf(' %g;',EactVect);
      pos.(sprintf('%s__Ea',vars{j})) = sprintf('[%s]',EactStr(1:end-1)); % delete ";" from end
    elseif Eact{j} ~= 0 && ~isnan(Eact{j})
      pos.(sprintf('%s__Ea',vars{j})) = sprintf('%g',1000*Eact{j});
    end
  end
  cellParams.function.pos = pos;
end

function [vars,vals,Eact] = parseData(data,fileName)
  % Prepare data by extracting only relevant part, deleting "Code Name"
  vars = data(:,2); vals = data(:,3); Eact = data(:,4);
  indVar = cellfun(@isstr,vars); % Index of strings
  indVar2 = logical(indVar.*indVar([end,1:end-1])); % Delete "Code Name"
  vars = vars(indVar2); vals = vals(indVar2); Eact = Eact(indVar2);
    
  indNum = cellfun(@isnumeric,vals);
  valsNum = sprintf('''@(x) (%g)''\n',[vals{indNum}]'); 
  vals(indNum) = eval(['{',valsNum,'}']);
  % Replace Inf, if any
  % strInf = strfind(vals,'Inf'); indInf = ~cellfun(@isempty,strInf);
  indInf = contains(vals,'Inf');
  vals(indInf) = {'@(x) (Inf)'};
  
  % Load vectors, if any
  vectInd = find(contains(vals,'['));
  if ~isempty(vectInd)
    for k = 1:length(vectInd)
      vectCell = vals(vectInd(k));
      vectCell = str2num(vectCell{:}); %#ok<ST2NM>
      vectCell = vectCell(:); % force to be column vector
      vectStr = sprintf(' %g;',vectCell);
      vals(vectInd(k)) = {['@(x) ([',vectStr(1:end-1),'])']};
    end
  end
  
  % Load tables, if any
  tableInd = find(contains(vals,'#'));
  if ~isempty(tableInd)
    for k = 1:length(tableInd)
      tabCell = vals(tableInd(k));
      tabCell = tabCell{1};
      tabName = tabCell(2:end);
      tableString = loadTable(fileName,tabName);
      vals(tableInd(k)) = {tableString};
    end
  end

  % Define all input arguments to be x and T
  for nf = 1:length(vals)
    fun_str = vals{nf};
    strFun_del1 = strfind(fun_str,'@(');
    strFun_del2 = strfind(fun_str,')');
    ind_del = find(strFun_del2>strFun_del1,1,'first');
    ind_arg = strFun_del1:strFun_del2(ind_del);
    arg_str = fun_str(ind_arg);
    vals{nf} = strrep(fun_str,arg_str,'@(x,T)');   
  end
end

% Combine Arrhenius relationship with reference parameters in function
function reg = addTemps(reg)
  R = 8.3144621;  % Universal gas constant
  f = fields(reg);
  for k = 1:length(f)
    fk = f{k};
    if length(fk)>4 && strcmp(fk(end-3:end),'__Ea')
      base = fk(1:end-4);
      if isfield(reg,base) % okay, we have a <param> and <param>__Ea combo
        baseExt = sprintf('.*exp(%s*(1/298.15-1/T)/%g)',reg.(fk),R);
        reg.(base) = [reg.(base) baseExt];
      end
      reg = rmfield(reg,fk);
    end
  end
end

% Converts from standard parameters to lumped parameters
% Output "cellData" is structure with STRINGS, not FUNCTIONS at this point
function [ cellData ] = convertStandardToLumped(cellDataStd)
  cellData = cellDataStd;
  cellData.standard = cellData.function;
  cellData = rmfield(cellData,'function'); % remove standard functions, replace
  cellDataFunc = cellStr2Fun(cellDataStd,0); % convert strings to functions
  cellDataFunc = cellDataFunc.function;

  % Reference temperature used in lumped-parameter computation
  % (only psi at this point).
  

  % Regions length
  LDL = cellDataFunc.DL.L();
  Lsep = cellDataFunc.sep.L();
  Lpos = cellDataFunc.pos.L(); 
  % Set up constants and calculations needed by transfer function
  F = cellData.const.F;
  R = cellData.const.R;
  TrefKelvin = cellData.const.T;
  t0plus  = cellDataFunc.const.t0plus();
  RsPos = cellDataFunc.pos.Rs();
  if isfield(cellDataFunc.pos,'Dsref')
    DsrefPos = cellDataFunc.pos.Dsref(); 
  else
    DsPos = cellDataFunc.pos.Ds();
  end
  sEpsPos  = cellDataFunc.pos.sEps();
  Acell = cellDataFunc.const.A();
  if exist('cellDataFunc.const.Rc()','var')
    Rc = cellDataFunc.const.Rc();
  else
    Rc = 0;
  end 
  asPos = 3*sEpsPos/RsPos;
  eEpsDL  = cellDataFunc.DL.eEps();
  eEpsSep  = cellDataFunc.sep.eEps();
  eEpsPos  = cellDataFunc.pos.eEps();
  De = cellDataFunc.const.De();
  % We could create "effective" De, but because we assume brugKappa and
  % brugDe are the same, the effective De is unused. Instead, we define psi
  % DeDL  = cellDataFunc.const.De() * eEpsDL^cellDataFunc.DL.brugDe();
  % DeSep  = cellDataFunc.const.De() * eEpsSep^cellDataFunc.sep.brugDe();
  % DePos  = cellDataFunc.const.De() * eEpsPos^cellDataFunc.pos.brugDe();
  ce0         = cellDataFunc.const.ce0();
  kappa       = cellDataFunc.const.kappa(ce0,0);
  kappaEffDL = kappa*(eEpsDL)^(cellDataFunc.DL.brugDeKappa());
  kappaEffSep = kappa*(eEpsSep)^(cellDataFunc.sep.brugDeKappa());
  kappaEffPos = kappa*(eEpsPos)^(cellDataFunc.pos.brugDeKappa());
  sigmaPos    = cellDataFunc.pos.sigma();
  sigmaEffPos = sigmaPos.*(sEpsPos).^cellDataFunc.pos.brugSigma();
  theta0Pos   = cellDataFunc.pos.theta0();
  theta100Pos = cellDataFunc.pos.theta100();
  csmaxPos    = cellDataFunc.pos.csmax();
  cLiNeg      = cellDataFunc.neg.cLi(); 
  gammaNeg    = cellDataFunc.neg.gamma();
  alphaNeg    = cellDataFunc.neg.alpha();
  alphaPos    = cellDataFunc.pos.alpha();
  knormNeg    = cellDataFunc.neg.knorm();
  knormPos    = cellDataFunc.pos.knorm();
  RfNeg       = cellDataFunc.neg.Rf();
  RfPos       = cellDataFunc.pos.Rf();
  dlnfdlnc    = cellDataFunc.const.dlnfdlnc();
  RdlNeg      = cellDataFunc.neg.Rdl();
  RdlPos      = cellDataFunc.pos.Rdl();
  nFPos       = cellDataFunc.pos.nF();
  nDLNeg      = cellDataFunc.neg.nDL();
  nDLPos      = cellDataFunc.pos.nDL();
  CdlNeg      = cellDataFunc.neg.Cdl();
  CdlPos      = cellDataFunc.pos.Cdl();
  wDLNeg      = cellDataFunc.neg.wDL();
  wDLPos      = cellDataFunc.pos.wDL();
  nESep       = cellDataFunc.sep.nE();
  nEDL        = cellDataFunc.DL.nE();

  % Generate Lumped parameters
  %package = @(x) sprintf('@(x,T)(%g)',x);
  package = @(x)['@(x,T) ([',sprintf(' %g;',x(:)),'])'];
  fn.neg.alpha     = package(alphaNeg);
  fn.neg.k0        = package(gammaNeg*Acell*knormNeg*F*ce0^(1-alphaNeg)*cLiNeg^(alphaNeg));
  fn.neg.Cdl       = package(CdlNeg*gammaNeg*Acell);
  fn.neg.Rdl       = package(RdlNeg/(gammaNeg*Acell));
  fn.neg.Rf        = package(RfNeg/(gammaNeg*Acell));
  fn.neg.nDL       = package(nDLNeg);
  fn.neg.wDL       = package(wDLNeg);
  fn.DL.kappa      = package(kappaEffDL*Acell/LDL);
  fn.DL.qe         = package(eEpsDL*ce0*Acell*LDL*F/(1-t0plus)/3600);
  fn.DL.nE         = package(nEDL);
  fn.sep.kappa     = package(kappaEffSep*Acell/Lsep);
  fn.sep.qe        = package(eEpsSep*ce0*Acell*Lsep*F/(1-t0plus)/3600);
  fn.sep.nE        = package(nESep);

  if isfield(cellData.standard.pos,'Uocp')
    % Non-MSMR model.
    fn.pos.Uocp      = cellData.standard.pos.Uocp;
    fn.pos.dUocp     = cellData.standard.pos.dUocp;
  else
    % MSMR model.
    fn.pos.X         = cellData.standard.pos.X;
    fn.pos.omega     = cellData.standard.pos.omega;
    fn.pos.U0        = cellData.standard.pos.U0;
  end
  fn.pos.sigma     = package(sigmaEffPos*Acell/Lpos);
  fn.pos.kappa     = package(kappaEffPos*Acell/Lpos);
  fn.pos.theta0    = package(theta0Pos);
  fn.pos.theta100  = package(theta100Pos);
  if isfield(cellDataFunc.pos,'Dsref') 
    fn.pos.Dsref   = package(DsrefPos/RsPos^2);
  else
    fn.pos.Ds      = package(DsPos/RsPos^2);
  end
  fn.pos.qe        = package(eEpsPos*ce0*Acell*Lpos*F/(1-t0plus)/3600);
  fn.pos.k0        = package(knormPos*asPos*Acell*Lpos);
  fn.pos.alpha     = package(alphaPos);
  fn.pos.Rf        = package(RfPos/(asPos*Acell*Lpos));
  fn.pos.nF        = package(nFPos);
  fn.pos.Cdl       = package(CdlPos*asPos*Acell*Lpos);
  fn.pos.nDL       = package(nDLPos);
  fn.pos.wDL       = package(wDLPos);
  fn.pos.Rdl       = package(RdlPos/(asPos*Acell*Lpos));
  fn.const.Q       = package(sEpsPos*Acell*Lpos*csmaxPos*...
                           abs(theta100Pos - theta0Pos)*F/3600);
  fn.const.Rc      = package(Rc);
  fn.const.psi     = package(F*De/kappa*ce0/(1-t0plus)/TrefKelvin);
  fn.const.kD      = package(2*R*(t0plus-1)/F*(1+dlnfdlnc));

  % Copy activation energies...
  reg = {'const','neg','DL','sep','pos'};
  for nr = 1:length(reg)
    regData = cellData.standard.(reg{nr});
    f = fields(regData);
    for k = 1:length(f)
      fk = f{k};
      if strcmp(fk,'kappa__Ea') % copy from "const" to other regions
        fn.neg.kappa__Ea = regData.kappa__Ea;
        fn.sep.kappa__Ea = regData.kappa__Ea;
        fn.pos.kappa__Ea = regData.kappa__Ea;
      end
      if length(fk)>4 && strcmp(fk(end-3:end),'__Ea')
        fnReg = fn.(reg{nr});
        if strcmp(fk(1:end-4),'knorm')
          fnReg.k0__Ea = regData.(fk);
        else
          fnReg.(fk) = regData.(fk);
        end
        fn.(reg{nr}) = fnReg;
      end
    end
  end

  % copy electrode SOC equation
  fn.pos.soc = cellData.standard.pos.soc;
  cellData.function = fn;
end

function tableStr = loadTable(fileName,tabName)
  tabData = readtable(fileName,'Sheet',tabName,'ReadVariableNames',false);
  tabData = table2array(tabData);
  tabData1 = mat2str(tabData(:,1));  % independent variable
  tabData2 = mat2str(tabData(:,2));  % dependent variable
  % Replace data by string function for table lookup
  tableStr = sprintf('@(x,T) (interp1(%s,%s,x,''linear'',''extrap''))',...
             tabData1,tabData2);
end

% Convert strings stored in the data structure to MATLAB inline functions
function cellData = cellStr2Fun(cellData,temps)
  reg = {'const','neg','DL','sep','pos'};
  for nr = 1:length(reg)
    if temps
      cellData.function.(reg{nr}) = addTemps(cellData.function.(reg{nr}));
    end
    y1 = cellData.function.(reg{nr});
    f1 = fields(y1);
    y2 = struct2cell(y1);
    fnInd = find(contains(y2,'@'));
    y3 = cellfun(@str2func,y2(fnInd),'UniformOutput',false);
    y4 = cell2struct(y3,f1(fnInd));
    cellData.function.(reg{nr}) = y4;
  end
end


% Make cell-level OCV function from electrode-level OCP functions
function cellParams = makeOCV(cellParams)
  SOC = 0:0.001:1;
  posTheta = cellParams.function.pos.soc(SOC,298.15);
  posOCP = cellParams.function.pos.Uocp(posTheta,298.15);
  cellOCV = posOCP ;
  SOCstr = sprintf(' %g;',SOC);
  SOCstr = SOCstr(1:end-1); % remove trailing ';'
  OCVstr = sprintf(' %g;',cellOCV);
  OCVstr = OCVstr(1:end-1); % ditto
  OCVfn = sprintf('@(x,T) interp1([%s],[%s],x,''linear'',''extrap'')',SOCstr,OCVstr);
  cellParams.function.const.Uocv = str2func(OCVfn);
end

% Make electrode-level OCP function from MSMR parameter values
function Uocp = makeOCP(reg)
  F = 96485.3365; % Faraday's constant
  R = 8.3144621;  % Universal gas constant
  SOCvec = 0.001:0.0001:0.999; % span entire electrode range    CHANGED 0.001

  U0       = reg.U0();
  X        = reg.X();
  omega    = reg.omega();
  
  % First, find relationship at Tref = 298.15K = 25 degC
  T = 298.15; f = F/(R*T);
  cost = @(U,Qdes) Qdes - sum(X./(1+exp(f*(U-U0)./omega)));
  OCVvec1 = 0*SOCvec;
  x0 = [0 5];
  for k = 1:length(OCVvec1)
    try
      OCVvec1(k) = fzero(cost,x0,[],SOCvec(k)); x0 = OCVvec1(k);
    catch
      OCVvec1(k) = NaN;
    end
  end
  
  % Next, find relationship at Tref = 299.15K = 26 degC
  T = 299.15; f = F/(R*T);
  cost = @(U,Qdes) Qdes - sum(X./(1+exp(f*(U-U0)./omega)));
  OCVvec2 = 0*SOCvec;
  x0 = [0 5];
  for k = 1:length(OCVvec2)
    try
      OCVvec2(k) = fzero(cost,x0,[],SOCvec(k)); x0 = OCVvec2(k);
    catch
      OCVvec2(k) = NaN;
    end
  end
  
  SOCstr = sprintf(' %g;',SOCvec);   SOCstr = ['[',SOCstr(1:end-1),']'];
  OCVstr1 = sprintf(' %g;',OCVvec1); OCVstr1 = ['[',OCVstr1(1:end-1),']'];
  OCVstr2 = sprintf(' %g;',OCVvec2-OCVvec1); 
  OCVstr2 = ['[',OCVstr2(1:end-1),']'];
  
  UocpFn = sprintf('@(x,T) interp1(%s,%s+(T(:)-298.15).*%s,x,''linear'',''extrap'')', ...
                   SOCstr,OCVstr1,OCVstr2);
  Uocp = str2func(UocpFn);         
end

% Make electrode-level dOCP function from MSMR parameter values
function dUocp = makedOCP(reg)
  F = 96485.3365; % Faraday's constant
  R = 8.3144621;  % Universal gas constant
  SOCvec = 0:0.001:1; % span entire electrode range
  dUdQ = 0*SOCvec;

  U0       = reg.U0();
  X        = reg.X();
  omega    = reg.omega();
  
  % Find dUdQ using Uocp at Tref
  T = 298.15; f = F/(R*T);
  Uocp = reg.Uocp(SOCvec,T); % get Uocp at Tref
  for k = 1:length(Uocp)
    x = X./(1+exp(f*(Uocp(k)-U0)./omega)); % get x_j for this U
    dQdU = -F/R*sum(x.*(1-x./X)./omega);
    dUdQ(k) = 1/dQdU;
  end
  
  SOCstr = sprintf(' %g;',SOCvec); SOCstr = ['[',SOCstr(1:end-1),']'];
  dOCVstr = sprintf(' %g;',dUdQ);  dOCVstr = ['[',dOCVstr(1:end-1),']'];
  
  % Make function, adding in temperature-dependence
  dUocpFn = sprintf('@(x,T) interp1(%s,T(:).*%s,x,''linear'',''extrap'')', ...
                    SOCstr,dOCVstr);
  dUocp = str2func(dUocpFn);    
end

% Make electrode-level d2OCP function from MSMR parameter values
function d2Uocp = maked2OCP(reg)
  F = 96485.3365; % Faraday's constant
  R = 8.3144621;  % Universal gas constant
  SOCvec = 0:0.001:1; % span entire electrode range
  d2U = 0*SOCvec;

  U0       = reg.U0();
  X        = reg.X();
  omega    = reg.omega();
  
  % Find dUdQ using Uocp at Tref
  T = 298.15; f = F/(R*T);
  Uocp = reg.Uocp(SOCvec,T); % get Uocp at Tref
  for k = 1:length(Uocp)
    g = exp(f*(Uocp(k)-U0)./omega); % g_j for this U
    x = X./(1+g); % x_j for this U
    dthetadU = -F/R*sum(x.*(1-x./X)./omega);
    d2U(k) = 1/dthetadU;
  end
  
  SOCstr = sprintf(' %g;',SOCvec); SOCstr = ['[',SOCstr(1:end-1),']'];
  dOCVstr = sprintf(' %g;',d2U);  dOCVstr = ['[',dOCVstr(1:end-1),']'];
  
  % Make function, adding in temperature-dependence
  dUocpFn = sprintf('@(x,T) interp1(%s,T(:).*%s,x,''linear'',''extrap'')', ...
                    SOCstr,dOCVstr);
  dUocp = str2func(dUocpFn);    
end
