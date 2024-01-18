% SIMFOM Simulate full-order model in COMSOL for specified iapp(t), T(t).
%
% This utility function executes a COMSOL model for specified 
% input-current and input-temperature profiles, starting at an initial 
% cell SOC, using the LiveLink for MATLAB interface. The returned COMSOL
% object can be saved to a file using mphsave.m, or loaded into the 
% COMSOL GUI using mphlaunch.m. 
%
% Usage 1: [FOM,FOMout] = simFOM(FOM,simData)
% Usage 2: [FOM,FOMout] = simFOM(...,'InputType',typestring)
% Usage 3: [FOM,FOMout] = simFOM(...,'DebugFlag',trueORfalse)
% Usage 4: [FOM,FOMout] = simFOM(...,'FixExchangeCurrent',true)
% Usage 5: [FOM,FOMout] = simFOM(...,'VcellOnly',true)
% 
% Inputs:
%   FOM = COMSOL object containing the full-order model, created by 
%     genFOM.m
%   simData = simulation profile loaded using "loadInput" or created
%     programmatically. A structure with the following required fields 
%     depending on the 'InputType' parameter:
%     - 'lut': SOC0, time, T, Iapp
%     - 'sin': SOC0, time, T, freq, mag
%     - 'pwm': SOC0, time, T, freq, mag, D, phase, squashTime
%     See below for a detailed description of each field.
%   InputType = type of input to use. typestring should be one of:
%     - 'lut' (default): regular lookup table / interpolation
%     - 'sin': sinusoidal input with specified amplitude and frequency
%     - 'pwm': pulse-width-modulated input
%   DebugFlag = logical value (true or false) indicating whether or not to
%     output status messages to the command window (default true)
%   FixExchangeCurrent = logical value (true or false) indicating whether
%     or not to fix the Butler-Volmer exchange current to initial values
%     (that is, exchange currents remain constant throughout the
%     simulation, remaining invariant with electrolyte and solid
%     concentrations). This option was added to emulate Murbach et. al. 
%     models of the second-harmonic respose where the exchange current
%     density is assumed to remain constant.
%   VcellOnly = when true, fetch only the cell voltage from the COMSOL
%     model, skipping all other variables (speeds up collecting results).
%
% Output:
%   FOM     = Updated COMSOL object containing the full-order model and all
%             simulated data  
%   FOMout  = Data structure with results from simulation
%
% simData structure fields:
%   SOC0 = the initial state-of-charge [%]
%   TSHIFT = amount by which to forward-shift input waveform in time [s]
%   time = a vector of simulation time values [s]
%   T = the scalar temperature or vector of temperature values
%     corresponding to the time vector [degC]
%   Iapp = a vector of applied current values corresponding to the
%     time vector [A]
%   freq = the cyclic frequency of the sinusoidal input or the PWM
%     waveform [Hz]
%   mag = the amplitude of the sinusoidal input or PWM waveform [A]
%   D = the duty ratio of the PWM input, 0.5 => 50% duty cycle [-]
%   phase = the phase-offset of the PWM input [rad]
%   squashTime = initial portion of applied current "squashed" [s]
%
% -- Note --
% Electrochemical Impedance Spectroscopy (EIS) Simulations:
%   EIS simulation executes much faster (and more accurately) in COMSOL 
%   when using the 'sin' input type rather than the default 'lut' 
%   (an analytic expression is employed instead of an interpolation 
%   function when using 'sin').
%
% -- Copyright --
% Copyright (c) 2021 by Gregory L. Plett and M. Scott Trimboli of the
% University of Colorado Colorado Springs (UCCS). 
% This work is licensed under a Creative Commons 
% Attribution-NonCommercial-ShareAlike 4.0 Intl. License, v. 1.0. It is 
% provided "as is", without express or implied warranty, for educational 
% and informational purposes only.
%
% -- Changelog --
% 2023.12.17 | Distinquish between Phis and PhisTilde=Phis-Vcell | Wes H.
% 2023.04.20 | Add option to retain only Vcell | Wesley Hileman
% 2023.04.09 | Modify for LMB | Wesley Hileman <whileman@uccs.edu>
% 2023.04.04 | Add 'InputType' parameter | Wesley Hileman <whileman@uccs.edu>

function [FOM,FOMout] = simFOM(genData,simData,varargin)
  import com.comsol.model.*      %#ok<NSTIMP>
  import com.comsol.model.util.* %#ok<NSTIMP>

  % Validate input parameters.
  iscomsolmodel = @(x)isscalar(x)&&strcmp(x.origin__,'genFOM');
  issimdata = @(x)isstruct(x);
  isinputtype = @(x)(isstring(x)||ischar(x))&&any(strcmpi(x,{'lut','sin','pwm'}));
  isdebugflag = @(x)isscalar(x)&&islogical(x);
  parser = inputParser;
  parser.addRequired('genData',iscomsolmodel);
  parser.addRequired('simData',issimdata);
  parser.addParameter('InputType','lut',isinputtype);
  parser.addParameter('DebugFlag',true,isdebugflag);
  parser.addParameter('FixExchangeCurrent',false,@(x)isscalar(x)&&islogical(x));
  parser.addParameter('VcellOnly',false,@(x)isscalar(x)&&islogical(x));
  parser.parse(genData,simData,varargin{:});
  p = parser.Results;  % structure of validated parameters
  debugFlag = p.DebugFlag;

  FOM = genData.FOM;
  cellModel = genData.cellModel;
  FOMout = [];
  
  % Shift time vector when COMSOL samples all output values by TSHIFT to
  % allow COMSOL to respond to abrupt changes in input current. 
  %
  % If TSHIFT=0, then COMSOL output does not have time to respond to step-
  % like instant changes in input current before sampling voltage (etc.), 
  % so the voltage change does not appear until the following time sample,
  % making the overall voltage profile appear delayed by one time sample.
  % That is, a nonzero TSHIFT is needed to capture feedthrough resistance
  % (ESR) effects. A TSHIFT of 0.1% of the sample period works; a shift of
  % 1% might be slightly better.
  if isfield(simData,'TSHIFT')
    TSHIFT = simData.TSHIFT;
  else
    % Default.
    TSHIFT = 0.0025;
  end
                 
  updateInputs;
  runStudy;
  collectResults;
  msg('\n');
  
  function updateInputs
    msg('Updating inputs in FOM...');
    tvec = simData.time(:);
    Tvec = simData.T(:);
    SOC0 = simData.SOC0;
    
    % Replace default input-current function
    FOM.func.remove('Ides'); % we'll replace the input function
    if strcmpi(p.InputType,'lut')
        % Lookup table.
        ivec = simData.Iapp(:);
        if length(tvec) == length(ivec),  ivec = ivec(1:end-1); end
        v1 = tvec(1:end-1); v2 = tvec(2:end); v3 = ivec;
        Ides = sprintf('''%g'' ''%g'' ''%g'';',[v1(:),v2(:),v3(:)]');
        Ides = eval(sprintf('{ %s }',Ides));
        FOM.func.create('Ides', 'Piecewise');
        FOM.func('Ides').label('Default current profile');
        FOM.func('Ides').set('funcname', 'inputCurrent');
        FOM.func('Ides').set('arg', 't');
        FOM.func('Ides').set('extrap', 'interior');
        FOM.func('Ides').set('smooth', 'contd2');
        FOM.func('Ides').set('smoothzone', '3E-7');
        FOM.func('Ides').set('pieces', Ides);
        FOM.func('Ides').set('argunit', 's');
        FOM.func('Ides').set('fununit', 'A');
    elseif strcmpi(p.InputType,'sin')
        % Sine wave.
        I = simData.mag;
        f0 = simData.freq;
        FOM.func.create('Ides', 'Analytic');
        FOM.func('Ides').label('Default current profile');
        FOM.func('Ides').set('funcname', 'inputCurrent');
        FOM.func('Ides').set('args', 't');
        FOM.func('Ides').set('fununit', 'A');
        FOM.func('Ides').setIndex('argunit', 's', 0);
        FOM.func('Ides').set('expr', sprintf('%g*sin(2*pi*%g*t)',I,f0));
    elseif strcmpi(p.InputType,'pwm')
        % Pulse-width modulated square-wave.
        tend = tvec(end);
        period = 1/simData.freq;
        D = simData.D;
        phase = simData.phase;
        mag = simData.mag;
        TS = simData.squashTime; % How much delay
        t = TS; VV = [0 TS 0];
        if phase ~= 0
          VV = [VV; TS TS+period/4 mag];
          t = TS + period/4;
        end
        while t <= tend
          VV = [VV; t t+period*D -mag; t+period*D t+period mag];
          t = t+period;
        end
        Ides = sprintf('''%g'' ''%g'' ''%g'';',VV');
        Ides = eval(sprintf('{ %s }',Ides));
        FOM.func.create('Ides', 'Piecewise');
        FOM.func('Ides').label('Default current profile');
        FOM.func('Ides').set('funcname', 'inputCurrent');
        FOM.func('Ides').set('arg', 't');
        FOM.func('Ides').set('extrap', 'interior');
        FOM.func('Ides').set('smooth', 'contd2');
        FOM.func('Ides').set('smoothzone', '3E-7');
        FOM.func('Ides').set('pieces', Ides);
        FOM.func('Ides').set('argunit', 's');
        FOM.func('Ides').set('fununit', 'A');
    end
    
    % Replace default input temperature function
    FOM.func.remove('Temperature');
    if length(Tvec)==1
        % Constant temperature.
        FOM.func.create('Temperature', 'Analytic');
        FOM.func('Temperature').set('funcname', 'inputTemperature');
        FOM.func('Temperature').set('args', 't');
        FOM.func('Temperature').set('argunit', 's');
        FOM.func('Temperature').set('fununit', 'K');
        FOM.func('Temperature').set('expr',sprintf('%g',Tvec+273.15));
    else
        % Time-variant temperature.
        if length(tvec) == length(Tvec),  Tvec = Tvec(1:end-1); end
        v1 = tvec(1:end-1); v2 = tvec(2:end); v3 = Tvec+273.15;
        Temperature = sprintf('''%g'' ''%g'' ''%g'';',[v1(:),v2(:),v3(:)]');
        Temperature = eval(sprintf('{ %s }',Temperature));
        FOM.func.create('Temperature', 'Piecewise');
        FOM.func('Temperature').set('funcname', 'inputTemperature');
        FOM.func('Temperature').set('arg', 't');
        FOM.func('Temperature').set('extrap', 'periodic');
        FOM.func('Temperature').set('smooth', 'contd2');
        FOM.func('Temperature').set('smoothzone', '3E-7');
        FOM.func('Temperature').set('pieces', Temperature);
        FOM.func('Temperature').set('argunit', 's');
        FOM.func('Temperature').set('fununit', 'K');
    end

    % Replace default initial SOC.
    z0 = SOC0/100;
    FOM.param.set('z0', num2str(z0), 'Initial cell SOC');
    if cellModel.MSMR
        theta0 = cellModel.function.pos.soc(z0);
        U0 = cellModel.function.pos.U0(0.5,Tvec(1)+273.15);
        X = cellModel.function.pos.X(0.5,Tvec(1)+273.15);
        omega = cellModel.function.pos.omega(0.5,Tvec(1)+273.15);
        f = cellModel.const.F/cellModel.const.R/(Tvec(1)+273.15);
        Uocp0  = fzero(@(U)(theta0-sum(X./(1+exp(f*(U-U0)./omega)))), mean(U0));
        FOM.physics('U').feature('init1').set('U', Uocp0);
        FOM.physics('phi_s').feature('init2').set('phi_s', Uocp0);
    end

    % Use fixed exchange current if requested by user.
    if p.FixExchangeCurrent
        % Fix the Li-Metal exchange current.
        FOM.variable('var1').set('i0','k0','Exchange current');

        % Fix the porous-electrode exchange current.
        if cellModel.MSMR
            % Fix the MSMR exchange currents i0j to constant initial values.
            x0 = X./(1+exp(f*(Uocp0-U0)./omega));
            for j = 1:length(U0)
                % Names of MSMR variables for use in constructing expressions.
                Xj = sprintf('X%d',j);
                omegaj = sprintf('omega%d',j);
                k0j = sprintf('k0%d',j);
                alphaj = sprintf('alpha%d',j);
                xjss = sprintf('%g',x0(j)); % Fix at initial value!
                i0j = sprintf('i0%d',j);
                FOM.variable('varsPos1d').set( ...
                    i0j,sprintf( ...
                        '%s*nicePow(%s,%s*%s)*nicePow(%s-%s,%s*(1-%s))', ...
                        k0j,xjss,omegaj,alphaj,Xj,xjss,omegaj,alphaj));
            end
        else
            % Fix exchange current i0 to constant initial value.
            FOM.variable('varsPos1d').set('i0', '(k0*nicePow(1-thetaInitPos,1-alpha)*nicePow(thetaInitPos,alpha))', 'Exchange current');
        end
    end
    
    % Replace default study values
    tvar = tvec+TSHIFT; 
    FOM.study('std1').feature('time').set('tlist', tvar);      
    FOM.sol('sol1').feature('t1').set('tlist', tvar);
  end   % updateInputs
  function runStudy 
    import com.comsol.model.* 
    import com.comsol.model.util.*    
    msg('\nRunning study...');
    if ~ismac, ModelUtil.showProgress(true); end    
    % mphlaunch(mphtags(FOM));    
    try
      FOM.sol('sol1').runAll;
    catch e
      warning(e.identifier,['COMSOL did not converge! Results computed ' ...
          'will be saved. Details: %s'],e.message);
    end    
  end      % runStudy
  function collectResults
    msg('\nCollecting results...');  

    if isfield(cellModel.function.const,'Rc')
        Rc = cellModel.function.const.Rc();
    else
        Rc = 0;
    end

    if p.VcellOnly
        % User requests the cell voltage only.
        input = mpheval(FOM,'Iapp');   
        FOMout.Iapp = input.d1(:,end);
        data_pos = mpheval(FOM,'phi_s','Edim',1,'Selection',3);
        FOMout.Vcell = data_pos.d1(:,end) - FOMout.Iapp*Rc;
        return;
    end

    % Negative electrode.
    vars = {'phise0','phi_e','ifn','idln','ifdln','eta'};
    data_neg = mpheval(FOM,vars,'selection',1,'edim','boundary');
    locs_neg0 = data_neg.p(:);
    negPhise0          = data_neg.d1;
    negIf0             = data_neg.d3;
    negIdl0            = data_neg.d4;
    negIfdl0           = data_neg.d5;
    negEta0            = data_neg.d6;

    % Electrolyte.
    data_elect = mpheval(FOM,{'theta_e','phi_e'},'Edim',1);
    locs_elect = data_elect.p(:);
    Phie   = data_elect.d2;
    Thetae = data_elect.d1;
    
    % Positive electrode.
    vars = {'phi_s','phi_e','theta_e','if','idl','thetass','ifdl','eta'};
    data_pos = mpheval(FOM,vars,'Edim',1,'Selection',3);
    locs_pos = data_pos.p(:);
    posPhis    = data_pos.d1;
    posPhise   = data_pos.d1 - data_pos.d2;
    posIf      = data_pos.d4;
    posIdl     = data_pos.d5;
    posThetass = data_pos.d6;
    posIfdl    = data_pos.d7;
    posEta     = data_pos.d8;

    % Store electrochemical variables.
    FOMout.Phie         = Phie;
    FOMout.PhieTilde    = Phie + negPhise0; % shift ground point to x=0+
    FOMout.Thetae       = Thetae;
    FOMout.Ifdl         = [negIfdl0 posIfdl];
    FOMout.If           = [negIf0 posIf];
    FOMout.Idl          = [negIdl0 posIdl];
    FOMout.Phis         = posPhis;
    FOMout.PhisTilde    = posPhis - data_pos.d1(:,end); % debias from cell cell voltage
    FOMout.Phise        = [negPhise0 posPhise];
    FOMout.Thetass      = posThetass;
    FOMout.Eta          = [negEta0 posEta];

    % Store x-locations.
    FOMout.xLocs.Phie       = locs_elect;
    FOMout.xLocs.PhieTilde  = locs_elect;
    FOMout.xLocs.Thetae     = locs_elect;
    FOMout.xLocs.Ifdl       = [locs_neg0; locs_pos];
    FOMout.xLocs.If         = [locs_neg0; locs_pos];
    FOMout.xLocs.Idl        = [locs_neg0; locs_pos];
    FOMout.xLocs.Phis       = locs_pos;
    FOMout.xLocs.PhisTilde  = locs_pos;
    FOMout.xLocs.Phise      = [locs_neg0; locs_pos];
    FOMout.xLocs.Thetass    = locs_pos;
    FOMout.xLocs.Eta        = [locs_neg0; locs_pos];
                          
    % Simulation control and primary output
    temp = mpheval(FOM,'T');
    FOMout.T = temp.d1(:,end);
    FOMout.TdegC = temp.d1(:,end)-273.15; % Degrees Celsius
    tvar = mpheval(FOM,'t');  tvar = tvar.d1(:,1);
    FOMout.time = tvar-TSHIFT;
    input = mpheval(FOM,'Ides');   
    FOMout.Ides = input.d1(:,end);
    input = mpheval(FOM,'Iapp');   
    FOMout.Iapp = input.d1(:,end);
    FOMout.Vcell = data_pos.d1(:,end) - FOMout.Iapp*Rc; 

    data_pos = mpheval(FOM,'thetasavg_pos');
    FOMout.posSOC = data_pos.d1;
  end % collectResults
  function msg(theText)
    if debugFlag
      fprintf(theText);
    end
  end % msg
end