function [etaTF,aux] = tfEta(s,locs,cellData)
%TFETA

cellData = tfCommon(s,cellData);
[tfVals, data] = cellData.common.tfData.h11.tfEta(locs);

etaTF = tfVals.';
aux.hfGain = data.hfGain.';
aux.dcGain = data.dcGain.';
aux.res0 = data.res0.';
aux.xLoc = locs(:).';
aux.names = strcat(data.regNames,'Phis').';

end