function [phiseTF,aux] = tfPhise(s,locs,cellData)
%TFPHISE

cellData = tfCommon(s,cellData);
[tfVals, data] = cellData.common.tfData.h11.tfPhiseStar(locs);

phiseTF = tfVals.';
aux.hfGain = data.hfGain.';
aux.dcGain = data.dcGain.';
aux.res0 = data.res0.';
aux.xLoc = locs(:).';
aux.names = strcat(data.regNames,'Phise').';

end
