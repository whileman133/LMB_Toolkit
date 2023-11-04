function [thetassTF,aux] = tfThetass(s,locs,cellData) 
%TFTHETASS

cellData = tfCommon(s,cellData);
[tfVals, data] = cellData.common.tfData.h11.tfThetassStar(locs);

thetassTF = tfVals.';
aux.hfGain = data.hfGain.';
aux.dcGain = data.dcGain.';
aux.res0 = data.res0.';
aux.xLoc = locs(:).';
aux.names = strcat(data.regNames,'Thetass').';
  
end  