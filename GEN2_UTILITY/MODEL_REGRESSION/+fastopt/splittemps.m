function modelVect = splittemps(model,modelspec)
%SPLITTEMPS Create model structures of multiplicty=1 at each temperature.

modelVect = fastopt.evaltemp(model,modelspec,modelspec.temps-273.15);

end