function [thetassTF,aux] = tfThetassInt(s,locs,cellData)
%TFTHETASSINT

cellData = tfCommon(s,cellData);
[tfVals, data] = cellData.common.tfData.h11.tfThetass(locs);

thetassTF = tfVals.';
aux.xLoc = locs(:).';
aux.names = strcat(data.regNames,'ThetassInt').';

end