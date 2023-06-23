% function [freqResp,hfGain,res0,tfData] = evalTF(tflist,s,cellData) 
% 
% Inputs:
%   tflist   = cell array contains transfer-function names with regions
%   s        = frequency samples to where TF is to be evaluated
%   cellData = structure with cell parameter constants
%
% Outputs:
%   freqResp = TFs frequency responses of each variable
%   hfGain   = TFs high frequency gains
%   res0     = residuals
%   tfData   = TFs "name" and "xLoc"
%
% This function evaluates transfer-functions at given "s" and "cellData".
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

function [freqResp,hfGain,res0,tfData] = evalTF(tflist,s,cellData) %#ok<INUSD>
  % Create storage for function outputs (cell array)
  freqResp = cell(length(tflist),1); % frequency response
  hfGain   = cell(length(tflist),1); % high-frequency gain
  res0     = cell(length(tflist),1); % res0
  names    = cell(length(tflist),1); % TF names
  xLoc     = cell(length(tflist),1); % TF x locations

  % Loop over all desired transfer-functions
  for ntf = 1:length(tflist)
    TFinfo    = strsplit(tflist{ntf},{'tf','(s,',',%s)'});
    var       = TFinfo{2}; % variable names
    xLoc{ntf} = str2num(TFinfo{end-1})'; %#ok<ST2NM> (DO NOT USE STR2DOUBLE)

    % Convert transfer-function name into m-file name that implements it
    tfName = sprintf('tf%s',lower(var)); tfName(3) = upper(tfName(3)); 
    if strcmpi(tfName(end-2:end),'int'), tfName(end-2) = 'I'; end

    % Execute the transfer function
    evalStr = sprintf('[freqResp{ntf}, aux] = %s(s,xLoc{ntf},cellData);',tfName);
    eval(evalStr);

    hfGain{ntf} = aux.hfGain;
    res0{ntf}   = aux.res0;
    names{ntf}  = aux.names; 
  end

  % Build correct structure to output (stack matrices)
  % row: locs; column: s
  freqResp = cell2mat(freqResp);
  hfGain   = cell2mat(hfGain);
  res0     = cell2mat(res0);
  tfData.names = cat(1,names{:});
  tfData.xLoc  = cell2mat(xLoc);
end