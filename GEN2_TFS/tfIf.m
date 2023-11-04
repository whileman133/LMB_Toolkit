function [ifTF,aux] = tfIf(s,locs,cellData)
%TFIF

cellData = tfCommon(s,cellData);
[tfVals, data] = cellData.common.tfData.h11.tfIf(locs);

ifTF = tfVals.';
aux.hfGain = data.hfGain.';
aux.dcGain = data.dcGain.';
aux.res0 = data.res0.';
aux.xLoc = locs(:).';
aux.names = strcat(data.regNames,'If').';

end

  