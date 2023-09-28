% function xraData = loadXRA(fileName)
% 
% Inputs:
%   fileName = name of Excel spreadsheet containing the XRA tuning values
%
% Outputs:
%   xraData  = a structure containing the information read from fileName,
%              organized such that genROM can create ROMS
%
% This function loads the XRA-control parameters from the Excel spreadsheet
% having filename "fileName".
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

function xraData = loadXRA(fileName)
  % Load Excel spreadsheet data
  [~,~,data]  = xlsread(fileName,'Parameters');

  % Scan file for markers delineating certain sections
  for ii = 1:length(data)
    switch lower( data{ii,1} )
      case {'environmental'}
          yEnv = {ii; 'env'};
      case {'xra method'}
          yDRA = {ii; 'xra'};
      case {'transfer function','tf','tf''s'}
          yTF = {ii;'tf'};
      otherwise
          % do nothing
    end
  end

  ind = [yEnv,yDRA,yTF,{length(data);0}];
  % Strip Off All NaN Cells
  for ii = 1:(length(ind)-1)
    for jj = (ind{1,ii}+2):(ind{1,ii+1})
      param{jj,ii} = data{jj,2}; %#ok<AGROW>
      value{jj,ii} = data{jj,3}; %#ok<AGROW>
      % unit{jj,ii}  = data{jj,4}; Not used right now, but maybe later
      if( isnan( param{jj,ii} ) ) 
        param{jj,ii} = [];  %#ok<AGROW>
      end
      if( isnan( value{jj,ii} ) ) 
        value{jj,ii} = [];  %#ok<AGROW>
      end
      % if( isnan( unit{jj,ii}  ) ) unit{jj,ii}  = []; end
    end
  end

  % Form the XRA Parameters
  for ii = 1:(length(ind)-1)
    for jj = 1:length(param)
      % For everything but the transfer functions
      if(~isempty( param{jj,ii} ) && (jj < yTF{1}) )
        if( ischar( value{jj,ii} ) )
          num = eval(value{jj,ii});
        else
          num = value{jj,ii};
        end
        xraData.(param{jj,ii}) = num;
      end

      % Form transfer function list
      if( strcmpi(param{jj,ii},'tflist') )
        kk = 1;
        tflist = {sprintf('%s',value{jj,ii})};
        while ( jj+kk <= size(value,1) && ~isempty( value{jj+kk,ii} ) )
          tflist = [tflist; 
                   {sprintf('%s',value{jj+kk,ii})} ]; %#ok<AGROW>
          kk = kk+1;
        end
      end
    end
  end
  xraData.tf = tflist;  

end