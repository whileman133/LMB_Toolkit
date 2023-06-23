% function FOM = addFOMvlim(FOM,Vmin,Vmax)
% 
% Inputs:
%   FOM     = COMSOL object containing the full-order model, created by
%             genFOM.m
%   Vmin    = Minimum discharge voltage, triggering CV mode
%   Vmax    = Maximum charge voltage, triggering CV mode
%
% Output:
%   FOM     = Updated COMSOL object containing the revised full-order model 
%
% This utility function modifies a COMSOL model to add a trigger that
% activates whenever the voltage exceeds Vmax; after that point in the
% simulation, the input current is modified to enforce constant-voltage
% output. This function uses the LiveLink for MATLAB interface. The
% returned COMSOL object can be saved to a file using mphsave.m, or loaded
% into the COMSOL GUI using mphlaunch.m. 
% 
% Copyright (c) 2021 by Gregory L. Plett and M. Scott Trimboli of the
% University of Colorado Colorado Springs (UCVmaxS). This work is licensed
% under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 Intl.
% License, v. 1.0. It is provided "as is", without express or implied
% warranty, for educational and informational purposes only.
%
% This file is provided as a supplement to: Plett, Gregory L. and Trimboli,
% M. Scott, "Battery Management Systems, Volume III, Physics-Based
% Methods," Artech House, 2021. 

function FOM = addFOMvlim(FOM,Vmin,Vmax)
  fprintf('Adding constant-voltage trigger(s) to FOM...\n');
  FOM.param.set('VsetMax', sprintf('%g[V]',Vmax),'Maximum permitted cell voltage');
  FOM.param.set('VsetMin', sprintf('%g[V]',Vmin),'Minimum permitted cell voltage');
  FOM.component('mod1d').cpl.create('aveop1', 'Average'); % need
  FOM.component('mod1d').cpl('aveop1').selection.geom('geom1d', 0);
  FOM.component('mod1d').cpl('aveop1').selection.set(4);

  FOM.component('mod1d').physics.create('Iapp', 'GlobalEquations', 'geom1d');
  FOM.component('mod1d').physics('Iapp').identifier('Iapp');
  FOM.component('mod1d').physics('Iapp').feature('ge1').set('DependentVariableQuantity', 'none');
  FOM.component('mod1d').physics('Iapp').feature('ge1').set('CustomDependentVariableUnit', 'A');
  FOM.component('mod1d').physics.create('CVmode', 'Events', 'geom1d');
  FOM.component('mod1d').physics('CVmode').identifier('CVmode');
  FOM.component('mod1d').physics('CVmode').create('ds1', 'DiscreteStates', -1);
  FOM.component('mod1d').physics('CVmode').create('is1', 'IndicatorStates', -1);
  FOM.component('mod1d').physics('CVmode').create('impl1', 'ImplicitEvent', -1);

  FOM.component('mod1d').physics('Iapp').label('Compute CVmode');
  FOM.component('mod1d').physics('Iapp').feature('ge1').set('name', 'Iapp');
  FOM.component('mod1d').physics('Iapp').feature('ge1').set('equation', ...
    'CVmax*CVmin*((Iapp-Ides)/1[A])+(1-CVmax)*(aveop1(phi_s)-Iapp*RcFN(1,T)-VsetMax)/VsetMax+(1-CVmin)*(aveop1(phi_s)-Iapp*RcFN(1,T)-VsetMin)/VsetMin');
  FOM.component('mod1d').physics('Iapp').feature('ge1').set('initialValueU', 'Ides');
  
  FOM.component('mod1d').physics('CVmode').feature('ds1').set('dim', 'CVmax');
  FOM.component('mod1d').physics('CVmode').feature('ds1').set('dimInit', 1);
  FOM.component('mod1d').physics('CVmode').feature('ds1').set('dimDescr', 'Exceeded Vmax flag');
  FOM.component('mod1d').physics('CVmode').feature('ds1').label('Initial States');
  FOM.component('mod1d').physics('CVmode').feature('is1').set('indDim', 'MaxV');
  FOM.component('mod1d').physics('CVmode').feature('is1').set('g', 'aveop1(phi_s)-Iapp*RcFN(1,T)-VsetMax');
  FOM.component('mod1d').physics('CVmode').feature('is1').set('dimInit', 0);
  FOM.component('mod1d').physics('CVmode').feature('is1').set('dimDescr', 'Exceeded Vmax?');
  FOM.component('mod1d').physics('CVmode').feature('is1').label('Condition');
  FOM.component('mod1d').physics('CVmode').feature('impl1').set('condition', 'MaxV>0');
  FOM.component('mod1d').physics('CVmode').feature('impl1').set('reInitName', 'CVmax');
  FOM.component('mod1d').physics('CVmode').feature('impl1').set('reInitValue', 0);
  FOM.component('mod1d').physics('CVmode').feature('impl1').label('Reinitialize CVmax');
  FOM.component('mod1d').physics('CVmode').feature('ds1').setIndex('dim', 'CVmin', 1, 0);
  FOM.component('mod1d').physics('CVmode').feature('ds1').setIndex('dimInit', 0, 1, 0);
  FOM.component('mod1d').physics('CVmode').feature('ds1').setIndex('dimDescr', '', 1, 0);
  FOM.component('mod1d').physics('CVmode').feature('ds1').setIndex('dimInit', 1, 1, 0);
  FOM.component('mod1d').physics('CVmode').feature('ds1').setIndex('dimDescr', 'Below Vmin flag', 1, 0);
  FOM.component('mod1d').physics('CVmode').feature('is1').setIndex('indDim', 'MinV', 1, 0);
  FOM.component('mod1d').physics('CVmode').feature('is1').setIndex('g', 0, 1, 0);
  FOM.component('mod1d').physics('CVmode').feature('is1').setIndex('dimInit', 0, 1, 0);
  FOM.component('mod1d').physics('CVmode').feature('is1').setIndex('dimDescr', '', 1, 0);
  FOM.component('mod1d').physics('CVmode').feature('is1').setIndex('g', 'aveop1(phi_s)-Iapp*RcFN(1,T)-VsetMin', 1, 0);
  FOM.component('mod1d').physics('CVmode').feature('is1').setIndex('dimDescr', 'Below Vmin?', 1, 0);
  FOM.component('mod1d').physics('CVmode').feature.duplicate('impl2', 'impl1');
  FOM.component('mod1d').physics('CVmode').feature('impl2').setIndex('reInitName', 'CVmin', 0, 0);
  FOM.component('mod1d').physics('CVmode').feature('impl2').set('condition', 'MinV<0');
  FOM.component('mod1d').physics('CVmode').feature('impl2').label('Reinitialize CVmin');