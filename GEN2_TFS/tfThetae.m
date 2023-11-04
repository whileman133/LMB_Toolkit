function [thetaeTF,aux] = tfThetae(s,locs,cellData)
%TFTHETAE

cellData = tfCommon(s,cellData);
[tfVals, data] = cellData.common.tfData.h11.tfThetae(locs);

thetaeTF = tfVals.';
aux.hfGain = data.hfGain.';
aux.dcGain = data.dcGain.';
aux.res0 = data.res0.';
aux.xLoc = locs(:).';
aux.names = strcat(data.regNames,'Thetae').';

end