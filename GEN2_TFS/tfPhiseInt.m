function [phiseTF,aux] = tfPhiseInt(s,locs,cellData)
%TFPHISEINT

cellData = tfCommon(s,cellData);
[tfVals, data] = cellData.common.tfData.h11.tfPhise(locs);

phiseTF = tfVals.';
aux.xLoc = locs(:).';
aux.names = strcat(data.regNames,'PhiseInt').';

end