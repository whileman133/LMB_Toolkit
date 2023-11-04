function [ifdlTF,aux] = tfIfdl(s,locs,cellData)
%TFIFDL

cellData = tfCommon(s,cellData);
[tfVals, data] = cellData.common.tfData.h11.tfIfdl(locs);

ifdlTF = tfVals.';
aux.hfGain = data.hfGain.';
aux.dcGain = data.dcGain.';
aux.res0 = data.res0.';
aux.xLoc = locs(:).';
aux.names = strcat(data.regNames,'Ifdl').';

end