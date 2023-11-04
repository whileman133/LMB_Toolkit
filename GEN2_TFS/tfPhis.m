function [phisTF,aux] = tfPhis(s,locs,cellData)
%TFPHIS

cellData = tfCommon(s,cellData);
[tfVals, data] = cellData.common.tfData.h11.tfPhisTilde(locs);

phisTF = tfVals.';
aux.hfGain = data.hfGain.';
aux.dcGain = data.dcGain.';
aux.res0 = data.res0.';
aux.xLoc = locs(:).';
aux.names = strcat(data.regNames,'Phis').';

end