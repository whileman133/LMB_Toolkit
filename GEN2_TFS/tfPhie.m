function [phieTF,aux] = tfPhie(s,locs,cellData) 
%TFPHIE

cellData = tfCommon(s,cellData);
[tfVals, data] = cellData.common.tfData.h11.tfPhieTilde(locs);

phieTF = tfVals.';
aux.hfGain = data.hfGain.';
aux.dcGain = data.dcGain.';
aux.res0 = data.res0.';
aux.xLoc = locs(:).';
aux.names = strcat(data.regNames,'Phie').';

end