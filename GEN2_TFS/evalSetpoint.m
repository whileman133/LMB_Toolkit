function cellModel = evalSetpoint(cellModel,s,soc,T)
%EVALSETPOINT

cellParams = getCellParams(cellModel,'TdegC',T-273.15,'socPct',soc*100);
cellParams.const.soc = soc;
cellParams.const.T = T;
fnames = fieldnames(cellParams);
for k = 1:length(fnames)
    fname = fnames{k};
    cellModel.(fname) = cellParams.(fname);
end

if ~isempty(s)
   cellModel.common = []; % force "tfCommon" to overwrite 
   cellModel = tfCommon(s,cellModel);
end

end % function