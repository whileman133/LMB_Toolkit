% function cellParams = evalSetpoint(cellParams,s,soc,T)
% 
% Inputs:
%   cellParams = data structure containing cell model parameter functions,
%                most likely loaded from an Excel spreadsheet
%   s          = vector of 1j*w where w is a vector of frequencies at which 
%                transfer functions will be evaluated. This argument can be
%                empty if these frequencies are unknown, but transfer-
%                function computations are sped up if a frequency vector is
%                specified at this time
%   soc        = cell state-of-charge (between 0 and 1) for the setpoint at
%                which the model is evaluated
%   T          = cell temperature (in K) for the setpoint at which the
%                model is evaluated
%
% Output:
%   cellParams = updated data structure containing original data as well as
%                specific values for parameters at the specified setpoint
%
% This utility function evaluates the cell-model parameter values for a
% specific state-of-charge and temperature setpoint. These parameter values
% are created by evaluating the functions in cellParams.function and are
% stored in cellParams.const, cellParams.neg, cellParams,sep, and
% cellParams.pos. If "s" is not empty, cellParams.common is created to
% store some computations common to evaluating all transfer functions -- if
% this operation is done once at this point it speeds up later transfer-
% function evaluations.
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

% Takes the cellParams.function.[ fields ] functions and evaluates at the
% specific given SOC and T setpoint to produce constant values 
% ("s" is a vector of frequencies to input to tfCommon)
function cellParams = evalSetpoint(cellParams,s,soc,T)
  F = fields(cellParams.function.const);
  for k = 1:length(F)
    fn = cellParams.function.const.(F{k});
    val = fn(soc,T);
    cellParams.const.(F{k}) = val;
  end

  F = fields(cellParams.function.neg);
  for k = 1:length(F)
    fn = cellParams.function.neg.(F{k});
    val = fn(soc,T);
    cellParams.neg.(F{k}) = val;
  end
  
  
  if isfield(cellParams.function,'dll')
      F = fields(cellParams.function.dll);
      SOC = cellParams.function.pos.soc(soc,T);
      for k = 1:length(F)
              fn = cellParams.function.dll.(F{k});
              val = fn(SOC,T);
              cellParams.dll.(F{k}) = val;
      end
      cellParams.pos.soc = SOC;
  elseif isfield(cellParams.function,'DL')
      F = fields(cellParams.function.DL);
      SOC = cellParams.function.pos.soc(soc,T);
      for k = 1:length(F)
              fn = cellParams.function.DL.(F{k});
              val = fn(SOC,T);
              cellParams.dll.(F{k}) = val;
      end
      cellParams.pos.soc = SOC;
  end
  
  F = fields(cellParams.function.sep);
  for k = 1:length(F)
    fn = cellParams.function.sep.(F{k});
    val = fn(soc,T);
    cellParams.sep.(F{k}) = val;
  end
  
  F = fields(cellParams.function.pos);
  posSOC = cellParams.function.pos.soc(soc,T);
  for k = 1:length(F)
    fn = cellParams.function.pos.(F{k});
    val = fn(posSOC,T);
    cellParams.pos.(F{k}) = val;
  end
  cellParams.pos.soc = posSOC;
  
  cellParams.const.soc = soc;
  cellParams.const.T = T;
  
  % SOC-dependent Ds

  if isfield(cellParams.pos,'Dsref')
    cellParams.pos.Ds = cellParams.pos.Dsref*cellParams.const.F/cellParams.const.R/T...
        *cellParams.pos.soc*(cellParams.pos.soc-1)*cellParams.pos.dUocp;
  end
  
  if ~isempty(s)
    cellParams.common = []; % 20200629 force "tfCommon" to overwrite 
    [C,L,J,Z,Rct] = tfCommon(s,cellParams);
    cellParams.common.s = s;
    cellParams.common.C = C;
    cellParams.common.L = L;
    cellParams.common.J = J;
    cellParams.common.Z = Z;
    cellParams.common.Rct = Rct;
    cellParams.common.ind = ['L1DL=1;L1s=2;L1p=3;L2p=4;' ...
      'c1DL=1;c2DL=2;c1s=3;c2s=4;c1p=5;c2p=6;c3p=7;c4p=8;'...
      'j1p=1;j2p=2;j3p=3;j4p=4;' ...
      'Zsep=1;Zsp=2;Isp=3;'];

  end
end