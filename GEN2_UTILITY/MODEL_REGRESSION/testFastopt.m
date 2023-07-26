% testFastopt.m

addpath(fullfile('..','..'));
TB.addpaths;

temps = [15 25];
params.x.a = fastopt.param('tempfcn','fix');  % +1
params.x.b = fastopt.param('tempfcn','Eact'); % +2
params.x.c = fastopt.param('tempfcn','lut');  % +2
params.y.d = fastopt.param('len',3,'tempfcn','fix');  % +3
params.y.e = fastopt.param('len',3,'tempfcn','Eact'); % +4
params.y.f = fastopt.param('len',3,'tempfcn','lut');  % +6
params.z.g = fastopt.param('fix',0);
spec = fastopt.modelspec(params,'temps',temps,'TrefdegC',25);
assert(spec.nvars==18,"Wrong number of free parameters!");

model1.x.a = 1;
model1.x.b = 2;
model1.x.bEact = 100;
model1.x.c = [3 4];
model1.y.d = [5; 6; 7];
model1.y.e = [8; 9; 10];
model1.y.eEact = 200;
model1.y.f = [
    11  14
    12  15
    13  16
];

vect1 = fastopt.pack(model1,spec);
model1Recovered = fastopt.unpack(vect1,spec);
vect1Recovered = fastopt.pack(model1Recovered,spec);
assert(all(vect1==vect1Recovered));

modelVect = fastopt.splittemps(model1Recovered,spec);