% GENFOM Create full-order COMSOL model from LMB cell model.
%
% This utility function creates a COMSOL object on the COMSOL server via
% the LiveLink for MATLAB interface, using the cell data structure loaded
% from an Excel spreadsheet using loadCellParams.m. The returned COMSOL
% object can be saved to a file using mphsave.m, or loaded into the COMSOL
% GUI using mphlaunch.m. Generally, however, the recommended workflow is to
% use simFOM.m as the next step to merge desired profiles of current and
% temperature versus time to the FOM object and to execute the simulation
% via the LiveLink.
%
% genData = genFOM(cellModel)
% 
% Inputs:
%   cellModel = cell-model data structure from LOADCELLPARAMS.m
%
% Output:
%   genData   = structure with the following fields:
%     name     : name of the model
%     FOM      : COMSOL object containing the full-order model.
%     cellModel: cell model used to generate the FOM
%
% -- Changelog --
% 2023.06.02 | nojac() for MSMR exchange currents, i0j | Wesley Hileman
% 2023.04.09 | Adapt to LMB, add support for MSMR Models | Wesley Hileman <whileman@uccs.edu>

function genData = genFOM(cellModel,varargin)
  import com.comsol.model.* 
  import com.comsol.model.util.*

  parser = inputParser;
  parser.addRequired('cellModel',@(x)isscalar(x)&&isstruct(x));
  parser.addParameter('DebugFlag',true,@(x)isscalar(x)&&islogical(x));
  parser.parse(cellModel,varargin{:});
  p = parser.Results; % strcuture of validated parameters
  debugFlag = p.DebugFlag;

  % Convert cell model to legacy Lumped-Parameter Model.
  % (genFOM still works with legacy models.)
  p.cellModel = convertCellModel(p.cellModel,'LLPM');

  FOM = [];
  regs = {'const','neg','dll','sep','pos'};

  % Convert U0, X, omega, k0, alpha to sets of scalar parameters for
  % each gallery (e.g., X1, X2, ...) and remove the vector parameters.
  % COMSOL can't handle vector parameters!
  % NOTE: if the model is not an MSMR model, no change occurs, but
  % paramnamesMSMR is still set for use in unit setting below.
  [cellModel, paramnamesMSMR] = explodeMSMR(p.cellModel);

  % Defining default values... these will be replaced by user profiles
  ivec = zeros(10,1); % rest for 10s
  tvec = (0:length(ivec)-1)';
  Tvec = 25*ones(length(ivec),1);
  initSOC = 0.5;

  createModel;
  setParameterValues;
  createFunctions;
  setVariableValues;
  createPhysics;
  createStudy;

  % normally do later, but for testing...
  % runStudy;
  
  createResults;
  organize;
  msg('\n');

  genData.name = 'Full-Order LMB Cell Model';
  genData.FOM = FOM;
  genData.cellModel = p.cellModel;
  genData.origin__ = 'genFOM';

  %% functions that do all the work
  function createModel
    import com.comsol.model.*
    import com.comsol.model.util.*
  
    msg('Creating model...');    
    try
      FOM = ModelUtil.create('LiMetalCell');
    catch
      error('Error creating FOM: COMSOL with MATLAB not running?')
    end
    FOM.label('Lithium-metal cell');
    FOM.author('Wesley A. Hileman');
    FOM.comments(sprintf(['Ideal lithium-metal cell dynamics (1d)\n\n' ...
      'This COMSOL FOM simulates a 1d DFN FOM of a LMB cell. ' ...
      'MSMR OCP and kinetics model is implemented.\n']));

    % Create 1d FOM
    FOM.modelNode.create('mod1d', false);
    FOM.modelNode('mod1d').label('1d cell FOM');
    FOM.modelNode('mod1d').identifier('mainDim');
    FOM.modelNode('mod1d').defineLocalCoord(false);
    FOM.modelNode('mod1d').curvedInterior(false);

    % Create 1d geometry
    FOM.geom.create('geom1d', 1);
    FOM.geom('geom1d').model('mod1d');
    FOM.mesh.create('mesh1d', 'geom1d');
    FOM.geom('geom1d').repairTolType('relative');
    FOM.geom('geom1d').feature.create('I1', 'Interval');
    FOM.geom('geom1d').feature('I1').set('intervals','many');
    FOM.geom('geom1d').feature('I1').set('p','range(0,1,3)');

    % Create 1d mesh
    FOM.mesh('mesh1d').feature.create('E1', 'Edge');
    % mesh for dead Li layer
    FOM.mesh('mesh1d').feature('E1').feature.create('D1', 'Distribution');
    FOM.mesh('mesh1d').feature('E1').feature('D1').label('Dead Li Layer');
    FOM.mesh('mesh1d').feature('E1').feature('D1').selection.set(1);
    FOM.mesh('mesh1d').feature('E1').feature('D1').set('type','predefined');
    FOM.mesh('mesh1d').feature('E1').feature('D1').set('elemcount','20');
    FOM.mesh('mesh1d').feature('E1').feature('D1').set('elemratio','0.5');
    FOM.mesh('mesh1d').feature('E1').feature('D1').set('method','geometric');
    % mesh for separator
    FOM.mesh('mesh1d').feature('E1').feature.create('D2', 'Distribution');
    FOM.mesh('mesh1d').feature('E1').feature('D2').label('Separator');
    FOM.mesh('mesh1d').feature('E1').feature('D2').selection.set(2);
    FOM.mesh('mesh1d').feature('E1').feature('D2').set('type','predefined');
    FOM.mesh('mesh1d').feature('E1').feature('D2').set('elemcount','20');
    FOM.mesh('mesh1d').feature('E1').feature('D2').set('elemratio','0.5');
    FOM.mesh('mesh1d').feature('E1').feature('D2').set('symmetric',true);
    FOM.mesh('mesh1d').feature('E1').feature('D2').set('reverse', true);
    % mesh for positive electrode
    FOM.mesh('mesh1d').feature('E1').feature.create('D3','Distribution');
    FOM.mesh('mesh1d').feature('E1').feature('D3').label('Positive electrode');
    FOM.mesh('mesh1d').feature('E1').feature('D3').selection.set(3);
    FOM.mesh('mesh1d').feature('E1').feature('D3').set('type','predefined');
    FOM.mesh('mesh1d').feature('E1').feature('D3').set('elemcount','20');
    FOM.mesh('mesh1d').feature('E1').feature('D3').set('elemratio','0.5');
    FOM.mesh('mesh1d').feature('E1').feature('D3').set('method','geometric');
    FOM.mesh('mesh1d').feature('E1').feature('D3').set('reverse', true);
    FOM.mesh('mesh1d').run;

    % Create 1d view
    FOM.component('mod1d').view('view1').axis.set('xmin', -0.15);
    FOM.component('mod1d').view('view1').axis.set('xmax',  3.15);
    FOM.component('mod1d').view('view1').label('View of 1d geometry');

    % Create 2d FOM
    FOM.modelNode.create('mod2d', false);
    FOM.modelNode('mod2d').label('2d electrode FOM');
    FOM.modelNode('mod2d').identifier('pseudoDim');
    FOM.modelNode('mod2d').defineLocalCoord(false);
    FOM.modelNode('mod2d').curvedInterior(false);

    % Create 2d geometry
    FOM.geom.create('geom2d', 2);
    FOM.geom('geom2d').label('2d electrode geometry');
    FOM.geom('geom2d').model('mod2d');
    FOM.mesh.create('mesh2d', 'geom2d');
    FOM.geom('geom2d').repairTolType('relative');
    FOM.geom('geom2d').create('Pos', 'Square');
    FOM.geom('geom2d').feature('Pos').name('Positive electrode');
    FOM.geom('geom2d').feature('Pos').set('pos', '2.0,0.0');
    FOM.geom('geom2d').feature('fin').set('repairtoltype', 'relative');

    % Create 2d mesh
    FOM.mesh('mesh2d').create('ftri1', 'FreeTri');
    FOM.mesh('mesh2d').feature('ftri1').create('dis1', 'Distribution');
    FOM.mesh('mesh2d').feature('ftri1').feature('dis1').selection.set([3]);
    FOM.mesh('mesh2d').label('Mesh 2');
    FOM.mesh('mesh2d').feature('size').set('hauto', 3);
    FOM.mesh('mesh2d').feature('size').set('custom', 'on');
    FOM.mesh('mesh2d').feature('ftri1').set('yscale', 0.5);
    FOM.mesh('mesh2d').feature('ftri1').feature('dis1').set('numelem', 600);
    FOM.mesh('mesh2d').run;

    % Create 2d view
    FOM.component('mod2d').view('view2').axis.set('xmin', -0.0625);
    FOM.component('mod2d').view('view2').axis.set('xmax',  3.0625);
    FOM.component('mod2d').view('view2').axis.set('ymin', -0.5);
    FOM.component('mod2d').view('view2').axis.set('ymax',  1.5);
    FOM.component('mod2d').view('view2').label('View of 2d geometry');

    FOM.variable.create('varsNeg1d');
    FOM.variable('varsNeg1d').label('Dead Li and Negative-Electrode variables');
    FOM.variable('varsNeg1d').model('mod1d');
    FOM.variable('varsNeg1d').selection.geom('geom1d', 1);
    FOM.variable('varsNeg1d').selection.set(1);

    FOM.variable.create('varsSep1d');
    FOM.variable('varsSep1d').label('Separator variables');
    FOM.variable('varsSep1d').model('mod1d');
    FOM.variable('varsSep1d').selection.geom('geom1d', 1);
    FOM.variable('varsSep1d').selection.set(2);

    FOM.variable.create('varsPos1d');
    FOM.variable('varsPos1d').label('Positive-electrode variables ');
    FOM.variable('varsPos1d').model('mod1d');
    FOM.variable('varsPos1d').selection.geom('geom1d', 1);
    FOM.variable('varsPos1d').selection.set(3);

    % Anode variables at x=0
    FOM.component('mod1d').variable.create('var1');
    FOM.component('mod1d').variable('var1').label('Anode interface');
    FOM.component('mod1d').variable('var1').selection.geom('geom1d', 0);
    FOM.component('mod1d').variable('var1').selection.set([1]);

    FOM.variable.create('varsPos2d');
    FOM.variable('varsPos2d').label('Positive-electrode variables');
    FOM.variable('varsPos2d').model('mod2d');
    FOM.variable('varsPos2d').selection.geom('geom2d', 2);
    FOM.variable('varsPos2d').selection.set([1]);
    
    FOM.cpl.create('linext2', 'LinearExtrusion', 'geom1d');
    FOM.cpl('linext2').selection.set(3);
    FOM.cpl('linext2').label('Linear extrusion, positive electrode');
    FOM.cpl('linext2').set('opname', 'pos_thetass');
    FOM.cpl('linext2').set('dstgeom', 'geom2d');
    FOM.cpl('linext2').set('dstframe', 'material');
    FOM.cpl('linext2').set('srcframe', 'material');
    FOM.cpl('linext2').selection('srcvertex1').set(3);
    FOM.cpl('linext2').selection('srcvertex2').set(4);
    FOM.cpl('linext2').selection('dstvertex1').set(2);
    FOM.cpl('linext2').selection('dstvertex2').set(4);

    % Extrusion to get value of porous-electrode variable at the surface of
    % the solid particle (r(tilde)=1).
    FOM.cpl.create('linext4', 'LinearExtrusion', 'geom2d');
    FOM.cpl('linext4').selection.geom('geom2d', 1);
    FOM.cpl('linext4').selection.set(3);
    FOM.cpl('linext4').label('Linear extrusion, positive electrode');
    FOM.cpl('linext4').set('opname', 'pos_ss');
    FOM.cpl('linext4').set('dstgeom', 'geom1d');
    FOM.cpl('linext4').set('dstframe', 'material');
    FOM.cpl('linext4').set('srcframe', 'material');
    FOM.cpl('linext4').selection('srcvertex1').set(2);
    FOM.cpl('linext4').selection('srcvertex2').set(4);
    FOM.cpl('linext4').selection('dstvertex1').set(3);
    FOM.cpl('linext4').selection('dstvertex2').set(4);

    FOM.cpl.create('linproj2', 'LinearProjection', 'geom2d');
    FOM.cpl('linproj2').selection.set(1);
    FOM.cpl('linproj2').selection('srcvertex1').set(1);
    FOM.cpl('linproj2').selection('srcvertex2').set(3);
    FOM.cpl('linproj2').selection('srcvertex3').set(2);
    FOM.cpl('linproj2').selection('dstvertex1').set(1);
    FOM.cpl('linproj2').selection('dstvertex2').set(3);

    FOM.cpl.create('genproj2', 'GeneralProjection', 'geom2d');
    FOM.cpl('genproj2').selection.set(1);
    FOM.cpl('genproj2').active(false);
    FOM.cpl.create('negint', 'Integration', 'geom1d');
    FOM.cpl('negint').selection.set(1);
    FOM.cpl('negint').label('Dead Li Layer integration');
    FOM.cpl.create('posint', 'Integration', 'geom1d');
    FOM.cpl('posint').selection.set(3);
    FOM.cpl('posint').label('Positive-electrode integration');

    % % for debugging purposes, display geometries in MATLAB plot windows
    % close all
    % figure; mphgeom(FOM,'geom1d')
    % figure; mphgeom(FOM,'geom2d')
    % 
    % % for dz0ebugging purposes, display meshes in MATLAB plot windows
    % figure; mphmesh(FOM,'mesh1d')
    % figure; mphmesh(FOM,'mesh2d')
    % return
  end
  function setParameterValues
    msg('\nSetting parameter values...');
    FOM.param.group.create('paramsCell');
    FOM.param('paramsCell').set('F', '96485.3365[C/mol]', 'Faraday''s constant');
    FOM.param('paramsCell').set('R', '8.3144621 [J/mol/K]', 'Gas constant');
    FOM.param('paramsCell').set('thetaInitPos', 'theta0PosFN(0.5,T)+z0*(theta100PosFN(0.5,T)-theta0PosFN(0.5,T))', 'Initial positive-electrode soc');
    FOM.param('paramsCell').label('Cell parameters');

    % default parameter group is "simulation control"
    FOM.param.set('z0', num2str(initSOC), 'Initial cell SOC');
    FOM.param.set('Ides', 'inputCurrent(t)', 'Applied current');
    FOM.param.set('Iapp', 'Ides', 'Limited current');
    FOM.param.set('T', 'inputTemperature(t)', 'Ambient temperature');
    FOM.param.set('xnorm', '1[m]');
    FOM.param.set('ynorm', '1[m]');
    FOM.param.set('Vnorm', '1[V]');
    FOM.param.set('molnorm', '1[mol]');
    FOM.param.set('Ahnorm', '1[A*h]');
    FOM.param.label('Simulation control'); % label doesn't "stick"
    
    for nr = 1:length(regs)
      reg_params = fields(cellModel.function.(regs{nr}));
      for np = 1:length(reg_params)
        genFn(cellModel,regs{nr},reg_params{np})
      end
    end
    
    function genFn(cellData,elec,param)
      fn = cellData.function.(elec).(param);
      fn_value = char(cellData.function.(elec).(param));
      % Do we need to generate table as well?
      [tab_ind1] = strfind(fn_value,'interp1(');
      if ~isempty(tab_ind1)
        obj = reverseFcn(fn);
        if strcmp(elec,'const')
          tab_name = sprintf('%sTAB',param);
        else
          tab_name = sprintf('%s%sTAB',param,camel(elec));
        end
        tab_value = [obj.value{1},obj.value{2}];

        FOM.func.create(tab_name, 'Interpolation');
        fTable = sprintf('''%g'' ''%g'';',tab_value.');
        fTable = eval(sprintf('{ %s }',fTable));
        FOM.func(tab_name).set('table', fTable);
        FOM.func(tab_name).set('extrap', 'linear');
        FOM.func(tab_name).set('funcname', tab_name);
        FOM.func(tab_name).set('argunit', '');
        FOM.func(tab_name).set('fununit','');
        % Replace interp1
        fn_value = sprintf('@(x,T)%s(x)',tab_name); 
      end
      % Create Function
      if strcmp(elec,'const')
        fn_name = sprintf('%sFN',param);
      else
        fn_name = sprintf('%s%sFN',param,camel(elec));
      end
      
      FOM.func.create(fn_name, 'Analytic');
      FOM.func(fn_name).name(fn_name);
      [fn_ind1] = strfind(fn_value,'@(');
      [fn_ind2] = strfind(fn_value,')');
      fn_ind2 = fn_ind2(find(fn_ind2>fn_ind1,1,'first'));
      % FOM.func(fn_name).set('expr', ['x/x*T/T*',fn_value(fn_ind2+1:end)]);
      theExpr = fn_value(fn_ind2+1:end);
      if strcmpi(theExpr,'(inf)')
        if strcmpi(param,'dUocp')
          % Okay since this controls solid diffusion equation
        else
          if ~strcmpi(param,'rdl') 
            warning('"Inf" found in parameter %s; setting to 1',param);
          end % don't warn for Rdl since Cdl==0 overpowers the "1"
          theExpr = '(1)';
        end
      end
      % Check to see if the expression contains ".*" or "./" or ".^"
      while true
        ind1 = strfind(theExpr,'.*'); if ~isempty(ind1), theExpr(ind1)=[]; end
        ind2 = strfind(theExpr,'./'); if ~isempty(ind2), theExpr(ind2)=[]; end
        ind3 = strfind(theExpr,'.^'); if ~isempty(ind3), theExpr(ind3)=[]; end
        if isempty(ind1) && isempty(ind2) && isempty(ind3), break, end
      end
      FOM.func(fn_name).set('expr',theExpr);
      FOM.func(fn_name).set('funcname', fn_name);

      if strcmp(elec,'const')
        FOM.func(fn_name).set('args', {'x' 'T'});
        FOM.func(fn_name).set('argunit', 'mol, K');
      else
        FOM.func(fn_name).set('args', {'x' 'T'});
        FOM.func(fn_name).set('argunit', 's/s,K');
      end
      
      switch param
        case {'sigma','kappa'}
          FOM.func(fn_name).set('fununit', 'S');
        case {'Uocp','dUocp','Uocv','Vmin','Vmax'}
          FOM.func(fn_name).set('fununit', 'V');
        case {'Rf','Rdl','Rc','Rct'}
          FOM.func(fn_name).set('fununit', 'V/A');
        case {'nF','nDL','nE','alpha','theta0','theta100','soc','Thetas'}
          FOM.func(fn_name).set('fununit', ''); 
        case {'Dsref','Ds'}
          FOM.func(fn_name).set('fununit', '1/s');
        case {'kD','psi'}  
          FOM.func(fn_name).set('fununit', 'V/K');
        case {'Q'}
          FOM.func(fn_name).set('fununit', 'A*h');
        case {'qe'}
          FOM.func(fn_name).set('fununit', 'A*h');
        case {'Cdl'}
          FOM.func(fn_name).set('fununit', 'F');
        case {'k0'}
          FOM.func(fn_name).set('fununit', 'A');
        case {'J'}
            FOM.func(fn_name).set('fununit', '');
        case {'R'}
            FOM.func(fn_name).set('fununit', 'J/mol/K');
        case {'F'}
            FOM.func(fn_name).set('fununit', 'C/mol');
        case {'T','Tref'}
            FOM.func(fn_name).set('fununit', 'K');
        case {'W'}
            FOM.func(fn_name).set('fununit', '');
        case {'tauW','tauF','tauDL'}
            FOM.func(fn_name).set('fununit', 's');
        case paramnamesMSMR.pos.U0
          FOM.func(fn_name).set('fununit','V');
        case paramnamesMSMR.pos.X
          FOM.func(fn_name).set('fununit','');
        case paramnamesMSMR.pos.omega
          FOM.func(fn_name).set('fununit','');
        case paramnamesMSMR.pos.k0
          FOM.func(fn_name).set('fununit','A');
        case paramnamesMSMR.pos.alpha
          FOM.func(fn_name).set('fununit','');
        otherwise
          warning('\n - Unknown units for param = %s',param); 
      end %switch
      function outStr = camel(inStr)
        outStr = lower(inStr);
        outStr(1) = upper(outStr(1));
      end
    end %genFn           
  end
  function createFunctions
    msg('\nCreating functions...');
    
    % (Default) input-current function
    if length(tvec) == length(ivec),  ivec = ivec(1:end-1); end
    v1 = tvec(1:end-1); v2 = tvec(2:end); v3 = ivec;
    Ides = sprintf('''%g'' ''%g'' ''%g'';',[v1(:),v2(:),v3(:)]');
    Ides = eval(sprintf('{ %s }',Ides));
    FOM.func.create('Ides', 'Piecewise');
    FOM.func('Ides').set('funcname', 'inputCurrent');
    FOM.func('Ides').set('arg', 't');
    FOM.func('Ides').set('extrap', 'periodic');
    FOM.func('Ides').set('smooth', 'contd2');
    FOM.func('Ides').set('smoothzone', '3E-7');
    FOM.func('Ides').set('pieces', Ides);
    FOM.func('Ides').set('argunit', 's');
    FOM.func('Ides').set('fununit', 'A');

    % (Default) input temperature function
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

    FOM.func.create('an1', 'Analytic');
    FOM.func('an1').label('Imaginary-free power');
    FOM.func('an1').set('funcname', 'nicePow');
    FOM.func('an1').set('expr', '(step(x)*abs(x))^(y)');    
    FOM.func('an1').set('args', {'x' 'y'});
    FOM.func('an1').set('plotargs', {'x' '-1' '1'; 'y' '0' '1'});

    FOM.func.create('step1', 'Step');
    FOM.func('step1').set('funcname', 'step');
    FOM.func('step1').set('smooth', 0.001); % 20200715    
  end
  function setVariableValues
    msg('\nSetting variable values...'); 

    % Lithium-metal electrode.
    FOM.variable('var1').set('PARAMETERS', '0', '--------------------------------------------');
    FOM.variable('var1').set('eta', '-phi_e-Rf*Iapp', 'Overpotential');
    FOM.variable('var1').set('SOLVE_if_idl', '0', '--------------------------------------------');
    FOM.variable('var1').set('ifn0', 'ifn - k0*(exp((1-alpha)*F*eta/(R*T)) - exp(-alpha*F*eta/(R*T)))');
    FOM.variable('var1').set('k0', 'k0NegFN(0.5,T)');
    FOM.variable('var1').set('alpha', 'alphaNegFN(0.5,T)');
    FOM.variable('var1').set('i0', '(k0*nicePow(theta_e,1-alpha))','Exchange current');
    FOM.variable('var1').set('Rf', 'RfNegFN(0.5,T)');
    FOM.variable('var1').set('Rdl', 'RdlNegFN(0.5,T)');
    FOM.variable('var1').set('Cdl', 'CdlNegFN(0.5,T)');
    FOM.variable('var1').set('idln0', 'idln - Cdl*(-phi_et-Rf*(ifnt+idlnt)) + Rdl*Cdl*idlnt');
    FOM.variable('var1').set('ifdln', 'ifn +idln');
    FOM.variable('var1').set('phise0', '-phi_e');

    % Dead-lithium layer.
    FOM.variable('varsNeg1d').set('PARAMETERS', '0', '--------------------------------------------');
    FOM.variable('varsNeg1d').set('kappa', 'kappaDllFN(theta_e,T)');
    FOM.variable('varsNeg1d').set('tauW', 'tauWDllFN(0.5,T)');
    FOM.variable('varsNeg1d').set('ifdl', '0[A]', 'No flux in dead Li layer');

    % Separator.
    FOM.variable('varsSep1d').set('PARAMETERS', '0', '--------------------------------------------');
    FOM.variable('varsSep1d').set('kappa', 'kappaSepFN(theta_e,T)');
    FOM.variable('varsSep1d').set('tauW', 'tauWSepFN(0.5,T)');
    FOM.variable('varsSep1d').set('ifdl', '0[A]', 'No flux in separator');

    % Porous electrode (pseudo dimension).
    FOM.variable('varsPos2d').set('theta0', 'theta0PosFN(0.5,T)');
    FOM.variable('varsPos2d').set('theta100', 'theta100PosFN(0.5,T)');
    FOM.variable('varsPos2d').set('soc', 'theta_s', 'Solid concentration ratio at any locations of the particle');
    FOM.variable('varsPos2d').set('socavg', '3*linproj2(theta_s*(y^2))/(ynorm^3)');
    FOM.variable('varsPos2d').set('thetaInit', 'thetaInitPos', 'Initial local state of charge');
    if isfield(cellModel.function.pos,'Dsref')
      if cellModel.MSMR
        FOM.variable('varsPos2d').set('Ds', 'nojac(-DsrefPosFN(0.5,T)*F/(R*T)*abs(soc*(soc-1))*dUdThetas)');
      else
        FOM.variable('varsPos2d').set('Ds', 'nojac(-DsrefPosFN(0.5,T)*F/(R*T)*abs(soc*(soc-1))*dUocpPosFN(soc,T))');
      end
    else
      FOM.variable('varsPos2d').set('Ds', 'DsPosFN(soc,T)');
    end

    % Porous electrode (main dimension).
    FOM.variable('varsPos1d').set('PARAMETERS', '0', '--------------------------------------------');
    FOM.variable('varsPos1d').set('sigma', 'sigmaPosFN(thetass,T)');
    FOM.variable('varsPos1d').set('kappa', 'kappaPosFN(theta_e,T)');
    FOM.variable('varsPos1d').set('tauW', 'tauWPosFN(0.5,T)');
    FOM.variable('varsPos1d').set('Rf', 'RfPosFN(thetass,T)');
    FOM.variable('varsPos1d').set('Rdl', 'RdlPosFN(thetass,T)');
    FOM.variable('varsPos1d').set('Cdl', 'CdlPosFN(thetass,T)');
    FOM.variable('varsPos1d').set('VARIABLES', '0', '--------------------------------------------');
    FOM.variable('varsPos1d').set('extrcpl_source_influx', '-if', 'Input lithium flux to pseudoDim');
    FOM.variable('varsPos1d').set('thetass', 'pseudoDim.pos_ss(extrcpl_source_thetass)', 'Surface concentration ratio');
    FOM.variable('varsPos1d').set('thetasavg', 'pseudoDim.pos_ss(socavg)', 'thetasavg(x,t)');
    FOM.variable('varsPos1d').set('thetasavg_pos', 'posint(thetasavg)/xnorm', 'thetasavg(t)');

    if cellModel.MSMR
        J = cellModel.function.pos.J();
        FOM.variable('varsPos1d').set('J','JposFN(thetass,T)','Number of MSMR galleries.');
        FOM.variable('varsPos1d').set('eta', 'phi_s-phi_e-Uss-Rf*ifdl', 'Surface Overpotential');
        FOM.variable('varsPos1d').set('SOLVE_if_idl', '0', '--------------------------------------------');
        FOM.variable('varsPos1d').set('if0', 'if-if_msmr', 'Faradaic rate (intercalation, MSMR)');
        FOM.variable('varsPos1d').set('idl0', 'idl-Cdl*(phi_st-phi_et-Rf*(ift+idlt))+Rdl*Cdl*idlt', 'Non-faradaic rate (double-layer)');
        FOM.variable('varsPos1d').set('ifdl', 'if+idl', 'Total interface lithium flux');
        FOM.variable('varsPos1d').set('SOLVE_MSMR', '0', '--------------------------------------------');
        FOM.variable('varsPos1d').set('Uss', 'pseudoDim.pos_ss(extrcpl_source_U)', 'Solid-surface potential.');
        for j = 1:J
            % Need xj in mainDim to compute the exchange current densities and
            % faradaic currents associated with each MSMR gallery!
            xj = sprintf('x%d',j);
            FOM.variable('varsPos1d').set( ...
                sprintf('%sss',xj), ...
                sprintf('pseudoDim.pos_ss(extrcpl_source_%s)',xj), ...
                sprintf('Gallery partial lithiation, j=%d',j));
        end

        % MSMR PARAMETERS - U0, X, Omega, k0, alpha
        namesparams = fieldnames(paramnamesMSMR.pos);
        for k = 1:length(namesparams)
            param = namesparams{k};
            names = paramnamesMSMR.pos.(param);
            for i = 1:length(names)
                name = names{i};
                % Make MSMR parameters available in both main dimension and the
                % pseudo dimension!
                FOM.variable('varsPos1d').set(name,sprintf('%sPosFN(0.5,T)',name));
                FOM.variable('varsPos2d').set(name,sprintf('%sPosFN(0.5,T)',name));
            end
        end
    
        % MSMR VARIABLES - gj(U), xj(U), i0j, ifj, theta_s, if_msmr=sum(ifj)
        exprThetas = '';
        exprDThetasDUocp = '';
        exprI0sum = '';
        exprIfsum = '';
        for j = 1:J
            % Names of MSMR variables for use in constructing expressions.
            U0j = sprintf('U0%d',j);
            Xj = sprintf('X%d',j);
            omegaj = sprintf('omega%d',j);
            k0j = sprintf('k0%d',j);
            alphaj = sprintf('alpha%d',j);
    
            % Function gj(U) for gallery j (throughout particle).
            gj = sprintf('g%d',j);
            FOM.variable('varsPos2d').set( ...
                gj,sprintf('exp((U-%s)*F/R/T/%s)',U0j,omegaj));
    
            % Partial lithiation of gallery j (throughout particle).
            xj = sprintf('x%d',j);
            FOM.variable('varsPos2d').set( ...
                xj,sprintf('%s/(1+%s)',Xj,gj));
    
            % Differential capacity (dxj/dU) for gallery j.
            dxjdU = sprintf('d%sdU',xj);
            FOM.variable('varsPos2d').set( ...
                dxjdU,sprintf('-(F/R/T)*%s*%s/%s/(1+%s)^2',Xj,gj,omegaj,gj));
    
            % Exchange current of gallery j (at particle surface).
            % ! use nojac() to avoid simulation errors in some cases.
            i0j = sprintf('i0%d',j);
            xjss = sprintf('x%dss',j);
            FOM.variable('varsPos1d').set( ...
                i0j,sprintf( ...
                    'nojac(%s*nicePow(%s,%s*%s)*nicePow(%s-%s,%s*(1-%s))*nicePow(theta_e,1-%s)/nicePow(%s/2,%s))', ...
                    k0j,xjss,omegaj,alphaj,Xj,xjss,omegaj,alphaj,alphaj,Xj,omegaj));
    
            % Faradiac current of gallery j (at particle surface).
            ifj = sprintf('if%d',j);
            FOM.variable('varsPos1d').set( ...
                ifj,sprintf( ...
                    '%s*( exp((1-%s)*F*eta/R/T) - exp(-%s*F*eta/R/T) )', ...
                    i0j,alphaj,alphaj));
    
            % Solid lithiation (sum of all xj).
            exprThetas = sprintf('%s + %s',exprThetas,xj);
    
            % Differential capacity (sum of all dxj/dU).
            exprDThetasDUocp = sprintf('%s + %s',exprDThetasDUocp,dxjdU);
    
            % Total exchange current (sum of all i0j).
            exprI0sum = sprintf('%s + %s',exprI0sum,i0j);
    
            % Total faradaic current (sum of all ifj).
            exprIfsum = sprintf('%s + %s',exprIfsum,ifj);
        end
        FOM.variable('varsPos1d').set('i0_sum',exprI0sum,'Total exchange current.');
        FOM.variable('varsPos1d').set('if_msmr',exprIfsum,'Total faradaic flux.');
        FOM.variable('varsPos2d').set('theta_s',exprThetas); % solid concentration
        FOM.variable('varsPos2d').set('dThetasdU',exprDThetasDUocp); % differential capacity
        FOM.variable('varsPos2d').set('dUdThetas','1/dThetasdU'); % reciprocal property holds for 1st derivative only (not 2nd derivative!)
    else % ~cellModel.MSMR
       FOM.variable('varsPos1d').set('J','0','Number of MSMR galleries.');
       FOM.variable('varsPos1d').set('k0', 'k0PosFN(thetass,T)');
       FOM.variable('varsPos1d').set('alpha', 'alphaPosFN(0.5,T)');
       FOM.variable('varsPos1d').set('i0', '(k0*nicePow(theta_e,1-alpha)*nicePow(1-thetass,1-alpha)*nicePow(thetass,alpha))', 'Exchange current');
       FOM.variable('varsPos1d').set('eta', 'phi_s-phi_e-UocpPosFN(thetass,T)-Rf*ifdl', 'Overpotential');
       FOM.variable('varsPos1d').set('SOLVE_if_idl', '0', '--------------------------------------------');
       FOM.variable('varsPos1d').set('if0', 'if-i0*(exp((1-alpha)*F*eta/(R*T))-exp(-alpha*F*eta/(R*T)))', 'Faradaic rate (intercalation)');
       FOM.variable('varsPos1d').set('idl0', 'idl-Cdl*(phi_st-phi_et-Rf*(ift+idlt))+Rdl*Cdl*idlt', 'Non-faradaic rate (double-layer)');
       FOM.variable('varsPos1d').set('ifdl', 'if+idl', 'Total interface lithium flux');
    end % if cellModel.MSMR

    FOM.variable.create('cplvar7');
    FOM.variable('cplvar7').label('Positive-electrode boundary');
    FOM.variable('cplvar7').model('mod2d');
    FOM.variable('cplvar7').selection.geom('geom2d', 1);
    FOM.variable('cplvar7').selection.set(3);
    FOM.variable('cplvar7').set('thetass_influx', 'mainDim.pos_thetass(extrcpl_source_influx)', 'Interfacial lithium flux at the particle surface (input from mainDim)');
    FOM.variable('cplvar7').set('extrcpl_source_thetass', 'theta_s', 'Output surface concentration back to mainDim');
    if cellModel.MSMR
        J = cellModel.function.pos.J();
        FOM.variable('cplvar7').set('extrcpl_source_U', 'U', 'Output surface potential back to mainDim');
        for j = 1:J
            xj = sprintf('x%d',j);
            FOM.variable('cplvar7').set( ...
                sprintf('extrcpl_source_%s',xj),xj, ...
                sprintf('Output partial lithiation %s back to mainDim',xj));
        end
    end
  end
  function createPhysics
    msg('\nCreating physics...');
    FOM.physics.create('phi_s', 'GeneralFormPDE', 'geom1d');
    FOM.physics('phi_s').model('mod1d');
    FOM.physics('phi_s').identifier('phi_s');
    FOM.physics('phi_s').field('dimensionless').field('phi_s');
    FOM.physics('phi_s').field('dimensionless').component({'phi_s'});
    FOM.physics('phi_s').prop('Units').set('DependentVariableQuantity', 'electricpotential');
    FOM.physics('phi_s').selection.set([3]);
    FOM.physics('phi_s').create('init2', 'init', 1);
    FOM.physics('phi_s').feature('init2').selection.set(3);
    FOM.physics('phi_s').create('cons1', 'Constraint', 0);
    FOM.physics('phi_s').feature('cons1').selection.set(1);
    FOM.physics('phi_s').create('flux1', 'FluxBoundary', 0);
    FOM.physics('phi_s').feature('flux1').selection.set(4);
    FOM.physics('phi_s').label('Solid charge conservation');
    FOM.physics('phi_s').prop('ShapeProperty').set('boundaryFlux', false);
    FOM.physics('phi_s').prop('Units').set('SourceTermQuantity', 'current');
    FOM.physics('phi_s').prop('Units').set('CustomSourceTermUnit', 'V');
    FOM.physics('phi_s').feature('gfeq1').set('f', '-ifdl');
    FOM.physics('phi_s').feature('gfeq1').set('Ga', '-sigma*phi_sx*(xnorm^2)');
    FOM.physics('phi_s').feature('gfeq1').set('da', 0);
    FOM.physics('phi_s').feature('gfeq1').label('General Form PDE');
    FOM.physics('phi_s').feature('init2').set('phi_s', 'UocpPosFN(thetaInitPos,T)');
    FOM.physics('phi_s').feature('cons1').set('R', 'phi_s');
    FOM.physics('phi_s').feature('flux1').set('g', '-Iapp*xnorm');

    FOM.physics.create('phi_e', 'GeneralFormPDE', 'geom1d');
    FOM.physics('phi_e').model('mod1d');
    FOM.physics('phi_e').identifier('phi_e');
    FOM.physics('phi_e').field('dimensionless').field('phi_e');
    FOM.physics('phi_e').field('dimensionless').component({'phi_e'});
    FOM.physics('phi_e').prop('Units').set('DependentVariableQuantity', 'electricpotential');
    FOM.physics('phi_e').label('Electrolyte charge conservation');
    FOM.physics('phi_e').prop('ShapeProperty').set('boundaryFlux', false);
    FOM.physics('phi_e').prop('Units').set('SourceTermQuantity', 'current');
    FOM.physics('phi_e').feature('gfeq1').set('f', 'ifdl');
    FOM.physics('phi_e').feature('gfeq1').set('Ga', '-kappa*(phi_ex-WFN(0.5,T)*psiFN(0.5,T)*T*theta_ex)*(xnorm^2)');
    FOM.physics('phi_e').feature('gfeq1').set('da', 0);
    FOM.physics('phi_e').feature('gfeq1').label('General Form PDE');
    FOM.physics('phi_e').feature('init1').set('phi_e', 0);
    FOM.physics('phi_e').create('flux1', 'FluxBoundary', 0);
    FOM.physics('phi_e').feature('flux1').selection.set([1]);
    FOM.physics('phi_e').feature('flux1').setIndex('g', '(ifn+idln)*xnorm', 0);

    % Implementing the OCP PDE for the 2D geommetry
    if cellModel.MSMR
        FOM.physics.create('U', 'ConvectionDiffusionEquation', 'geom2d');
        FOM.physics('U').model('mod2d');
        FOM.physics('U').identifier('U');
        FOM.physics('U').field('dimensionless').field('U');
        FOM.physics('U').prop('Units').set('DependentVariableQuantity', 'none');
        FOM.physics('U').prop('Units').set('CustomDependentVariableUnit', 'V');
        FOM.physics('U').create('flux1', 'FluxBoundary', 1);
        FOM.physics('U').feature('flux1').selection.set([3]);
        FOM.physics('U').label('Solid diffusion');
        FOM.physics('U').prop('ShapeProperty').set('boundaryFlux', false);
        FOM.physics('U').prop('Units').set('CustomSourceTermUnit', '1/s');
        FOM.physics('U').feature('cdeq1').set('f', 0);
        FOM.physics('U').feature('cdeq1').set('c', {'0' '0' '0' '(y^2)*DsrefPosFN(0.5,T)*(F/R/T)*soc*(soc-1)'});
        FOM.physics('U').feature('cdeq1').set('da', '((y/ynorm)^2)*dThetasdU');
        FOM.physics('U').feature('cdeq1').label('Solid diffusion');
        FOM.physics('U').feature('init1').set('U', '0');
        FOM.physics('U').feature('flux1').set('g', 'thetass_influx/(10800[s/h]*QFN(0.5,T))*abs(theta100-theta0)*ynorm'); 
    else
        FOM.physics.create('theta_s', 'ConvectionDiffusionEquation', 'geom2d');
        FOM.physics('theta_s').model('mod2d');
        FOM.physics('theta_s').identifier('theta_s');
        FOM.physics('theta_s').field('dimensionless').field('theta_s');
        FOM.physics('theta_s').prop('Units').set('DependentVariableQuantity', 'none');
        FOM.physics('theta_s').prop('Units').set('CustomDependentVariableUnit', 1);
        FOM.physics('theta_s').create('flux1', 'FluxBoundary', 1);
        FOM.physics('theta_s').feature('flux1').selection.set([3]);
        FOM.physics('theta_s').label('Solid diffusion');
        FOM.physics('theta_s').prop('ShapeProperty').set('boundaryFlux', false);
        FOM.physics('theta_s').prop('Units').set('CustomSourceTermUnit', '1/s');
        FOM.physics('theta_s').feature('cdeq1').set('f', 0);
        FOM.physics('theta_s').feature('cdeq1').set('c', {'0' '0' '0' '(y^2)*Ds'});
        FOM.physics('theta_s').feature('cdeq1').set('da', '((y/ynorm)^2)');
        FOM.physics('theta_s').feature('cdeq1').label('Solid diffusion');
        FOM.physics('theta_s').feature('init1').set('theta_s', 'thetaInit');
        FOM.physics('theta_s').feature('flux1').set('g', 'thetass_influx/(10800[s/h]*QFN(0.5,T))*abs(theta100-theta0)*ynorm'); 
    end

    FOM.physics.create('theta_e', 'ConvectionDiffusionEquation', 'geom1d');
    FOM.physics('theta_e').model('mod1d');
    FOM.physics('theta_e').identifier('theta_e');
    FOM.physics('theta_e').field('dimensionless').field('theta_e');
    FOM.physics('theta_e').prop('Units').set('DependentVariableQuantity', 'none');
    FOM.physics('theta_e').prop('Units').set('CustomDependentVariableUnit', 1);
    FOM.physics('theta_e').create('zflx2', 'ZeroFluxBoundary', 0);
    FOM.physics('theta_e').feature('zflx2').selection.set([4]);
    FOM.physics('theta_e').create('flux1', 'FluxBoundary', 0);
    FOM.physics('theta_e').feature('flux1').setIndex('g', 'Iapp*xnorm', 0);
    FOM.physics('theta_e').feature('flux1').selection.set([1]);
    FOM.physics('theta_e').label('Electrolyte diffusion');
    FOM.physics('theta_e').prop('ShapeProperty').set('boundaryFlux', false);
    FOM.physics('theta_e').prop('Units').set('CustomSourceTermUnit', 'A');
    FOM.physics('theta_e').feature('cdeq1').set('f', 'ifdl');
    FOM.physics('theta_e').feature('cdeq1').set('c', 'psiFN(0.5,T)*T*kappa*(xnorm^2)');
    FOM.physics('theta_e').feature('cdeq1').set('da', 'tauW*kappa*psiFN(0.5,T)*T');
    FOM.physics('theta_e').feature('cdeq1').label('Convection-Diffusion Equation');
    FOM.physics('theta_e').feature('init1').set('theta_e', 1);

    FOM.physics.create('if', 'DomainODE', 'geom1d');
    FOM.physics('if').model('mod1d');
    FOM.physics('if').identifier('if');
    FOM.physics('if').field('dimensionless').field('if');
    FOM.physics('if').field('dimensionless').component({'if'});
    FOM.physics('if').prop('Units').set('DependentVariableQuantity', 'none');
    FOM.physics('if').prop('Units').set('CustomDependentVariableUnit', 'A');
    FOM.physics('if').selection.set([3]);
    FOM.physics('if').label('Faradaic interface lithium flux');
    FOM.physics('if').prop('Units').set('CustomSourceTermUnit', 'A');
    FOM.physics('if').feature('dode1').set('f', 'if0');
    FOM.physics('if').feature('dode1').set('da', 0);

    FOM.physics.create('idl', 'DomainODE', 'geom1d');
    FOM.physics('idl').model('mod1d');
    FOM.physics('idl').identifier('idl');
    FOM.physics('idl').field('dimensionless').field('idl');
    FOM.physics('idl').field('dimensionless').component({'idl'});
    FOM.physics('idl').prop('Units').set('DependentVariableQuantity', 'none');
    FOM.physics('idl').prop('Units').set('CustomDependentVariableUnit', 'A');
    FOM.physics('idl').selection.set([3]);
    FOM.physics('idl').label('Non-faradaic molar flow rate');
    FOM.physics('idl').prop('Units').set('CustomSourceTermUnit', 'A');
    FOM.physics('idl').feature('dode1').set('f', 'idl0');
    FOM.physics('idl').feature('dode1').set('da', 0);


    FOM.physics.create('ifn', 'BoundaryODE', 'geom1d');
    FOM.physics('ifn').model('mod1d');
    FOM.physics('ifn').identifier('ifn');
    FOM.physics('ifn').label('Faradaic interface lithium flux at anode');
    FOM.physics('ifn').prop('Units').set('DependentVariableQuantity', 'current');
    FOM.physics('ifn').prop('Units').set('SourceTermQuantity', 'current');
    FOM.physics('ifn').feature('dode1').set('f', 'ifn0');
    FOM.physics('ifn').feature('dode1').set('da', 0);
    FOM.physics('ifn').selection.set([1]);
    FOM.physics('ifn').field('dimensionless').component(1, 'ifn');



    FOM.physics.create('idln', 'BoundaryODE', 'geom1d');
    FOM.physics('idln').model('mod1d');
    FOM.physics('idln').identifier('idln');
    FOM.physics('idln').label('Non-faradaic molar flow rate at anode');
    FOM.physics('idln').selection.set([1]);
    FOM.physics('idln').prop('Units').set('DependentVariableQuantity', 'current');
    FOM.physics('idln').prop('Units').set('SourceTermQuantity', 'current');
    FOM.physics('idln').field('dimensionless').component(1, 'idln');
    FOM.physics('idln').feature('dode1').setIndex('da', 0, 0);
    FOM.physics('idln').feature('dode1').setIndex('f', 'idln0', 0);


  end
  function createStudy
    msg('\nCreating study...');
    FOM.study.create('std1');
    FOM.study('std1').create('time', 'Transient');
    FOM.study('std1').label('Study');
    FOM.study('std1').feature('time').set('tlist', 'range(0.01,1,9.01)');
    FOM.study('std1').feature('time').set('usertol', true);
    FOM.study('std1').feature('time').set('rtol', '0.0001');   
    
    FOM.sol.create('sol1');
    FOM.sol('sol1').study('std1');
    FOM.sol('sol1').attach('std1');
    FOM.sol('sol1').create('st1', 'StudyStep');
    FOM.sol('sol1').create('v1', 'Variables');
    FOM.sol('sol1').create('t1', 'Time');
    FOM.sol('sol1').feature('t1').create('fc1', 'FullyCoupled');
    FOM.sol('sol1').feature('t1').create('d1', 'Direct');
    FOM.sol('sol1').label('Solver');
    FOM.sol('sol1').feature('v1').set('control', 'user');
    FOM.sol('sol1').feature('v1').set('resscalemethod', 'auto');
    FOM.sol('sol1').feature('v1').set('clist', {'range(0,1,9)' '1.5[s]'});
    FOM.sol('sol1').feature('v1').feature('mainDim_phi_e').label('mainDim.phi_e');
    FOM.sol('sol1').feature('v1').feature('mainDim_phi_s').label('mainDim.phi_s');
    FOM.sol('sol1').feature('t1').set('tlist', 'range(0,1,9)');
    FOM.sol('sol1').feature('t1').set('rtol', '0.0001');  
    FOM.sol('sol1').feature('t1').set('atolglobalmethod', 'unscaled');
    FOM.sol('sol1').feature('t1').set('atolglobalvaluemethod', 'manual');
    FOM.sol('sol1').feature('t1').set('atolglobal', '0.0001');  
    FOM.sol('sol1').feature('t1').set('ewtrescale', false);
    if cellModel.MSMR
        FOM.sol('sol1').feature('t1').set('atolvaluemethod', ...
          {'mainDim_phi_e' 'manual' 'mainDim_phi_s' 'manual' ...
          'pseudoDim_U' 'manual' 'mainDim_theta_e' 'manual' ...
          'mainDim_if' 'factor' 'mainDim_idl' 'factor'});
    else
        FOM.sol('sol1').feature('t1').set('atolvaluemethod', ...
          {'mainDim_phi_e' 'manual' 'mainDim_phi_s' 'manual' ...
          'pseudoDim_theta_s' 'manual' 'mainDim_theta_e' 'manual' ...
          'mainDim_if' 'factor' 'mainDim_idl' 'factor'});
    end
    FOM.sol('sol1').feature('t1').set('tstepsbdf', 'strict');
    FOM.sol('sol1').feature('t1').set('initialstepbdf', '0.0010');
    FOM.sol('sol1').feature('t1').set('maxstepconstraintbdf', 'const');
    FOM.sol('sol1').feature('t1').set('bwinitstepfrac', '1.0');
    FOM.sol('sol1').feature('t1').set('storeudot', false);
    FOM.sol('sol1').feature('t1').feature('dDef').set('ooc', false);
    FOM.sol('sol1').feature('t1').feature('dDef').set('rhob', 400);
    FOM.sol('sol1').feature('t1').feature('fcDef').set('probesel', 'manual');
    FOM.sol('sol1').feature('t1').feature('fc1').active(false);
    FOM.sol('sol1').feature('t1').feature('fc1').set('damp', '1.0');
    FOM.sol('sol1').feature('t1').feature('fc1').set('ratelimitactive', true);
    FOM.sol('sol1').feature('t1').feature('fc1').set('probesel', 'manual');
    FOM.sol('sol1').feature('t1').feature('d1').set('linsolver', 'pardiso');
    FOM.sol('sol1').feature('t1').feature('d1').set('ooc', false);
    FOM.sol('sol1').feature('t1').feature('d1').set('pardreorder', 'ndmt');
    FOM.sol('sol1').feature('t1').feature('d1').set('pardmtsolve', false);
    FOM.sol('sol1').feature('t1').feature('d1').set('errorchk', false);
  end
  function runStudy
    msg('\nRunning study...');
    FOM.sol('sol1').runAll;
  end
  function createResults
    msg('\nCreating results...');
    FOM.result.dataset('dset1').label('1d solution pulses');
    FOM.result.dataset('dset2').label('2d solution pulses');

    FOM.result.create('pg3', 'PlotGroup1D');
    FOM.result('pg3').create('lngr1', 'LineGraph');
    FOM.result('pg3').feature('lngr1').set('data', 'dset1');
    FOM.result('pg3').feature('lngr1').set('xdata', 'expr');
    FOM.result('pg3').feature('lngr1').selection.all;
    FOM.result('pg3').label('phi_s vs x');
    FOM.result('pg3').set('titletype', 'manual');
    FOM.result('pg3').set('title', 'Solid potential <B> \phi<sub>s</sub></B>');
    FOM.result('pg3').set('xlabel', 'Normalized x-coordinate (unitless)');
    FOM.result('pg3').set('xlabelactive', true);
    FOM.result('pg3').set('ylabel', 'Potential (V)');
    FOM.result('pg3').set('ylabelactive', true);
    FOM.result('pg3').feature('lngr1').set('xdataexpr', 'x');
    FOM.result('pg3').feature('lngr1').set('xdataunit', 'm');
    FOM.result('pg3').feature('lngr1').set('xdatadescr', 'x-coordinate');
    FOM.result('pg3').feature('lngr1').set('smooth', 'none');
    FOM.result('pg3').feature('lngr1').set('resolution', 'normal');
    
    FOM.result.create('pg4', 'PlotGroup1D');
    FOM.result('pg4').create('lngr1', 'LineGraph');
    FOM.result('pg4').feature('lngr1').set('xdata', 'expr');
    FOM.result('pg4').feature('lngr1').selection.all;
    FOM.result('pg4').feature('lngr1').set('expr', 'phi_e');
    FOM.result('pg4').label('phi_e vs x');
    FOM.result('pg4').set('titletype', 'manual');
    FOM.result('pg4').set('title', 'Electrolyte potential <B>\phi<sub>e</sub></B>');
    FOM.result('pg4').set('xlabel', 'Normalized x-coordinate (unitless)');
    FOM.result('pg4').set('xlabelactive', true);
    FOM.result('pg4').set('ylabel', 'Potential (V)');
    FOM.result('pg4').set('ylabelactive', true);
    FOM.result('pg4').feature('lngr1').set('xdataexpr', 'x');
    FOM.result('pg4').feature('lngr1').set('xdataunit', 'm');
    FOM.result('pg4').feature('lngr1').set('xdatadescr', 'x-coordinate');
    FOM.result('pg4').feature('lngr1').set('smooth', 'none');
    FOM.result('pg4').feature('lngr1').set('resolution', 'normal');

    FOM.result.create('pg5', 'PlotGroup1D');
    FOM.result('pg5').create('lngr1', 'LineGraph');
    FOM.result('pg5').feature('lngr1').set('xdata', 'expr');
    FOM.result('pg5').feature('lngr1').selection.all;
    FOM.result('pg5').feature('lngr1').set('expr', 'theta_e');
    FOM.result('pg5').label('theta_e vs x');
    FOM.result('pg5').set('titletype', 'manual');
    FOM.result('pg5').set('title', 'Electrolyte concentration ratio<B> \theta<sub>e</sub></B>');
    FOM.result('pg5').set('xlabel', 'Normalized x-coordinate (unitless)');
    FOM.result('pg5').set('xlabelactive', true);
    FOM.result('pg5').set('ylabel', 'Concentration ratio (unitless)');
    FOM.result('pg5').set('ylabelactive', true);
    FOM.result('pg5').feature('lngr1').set('xdataexpr', 'x');
    FOM.result('pg5').feature('lngr1').set('xdataunit', 'm');
    FOM.result('pg5').feature('lngr1').set('xdatadescr', 'x-coordinate');
    FOM.result('pg5').feature('lngr1').set('smooth', 'none');
    FOM.result('pg5').feature('lngr1').set('resolution', 'normal');

    FOM.result.create('pg6', 'PlotGroup2D');
    FOM.result('pg6').create('surf1', 'Surface');
    FOM.result('pg6').label('theta_s vs (x,r)');
    FOM.result('pg6').set('titletype', 'custom');
    FOM.result('pg6').set('typeintitle', false);
    FOM.result('pg6').set('descriptionintitle', false);
    FOM.result('pg6').set('suffixintitle', ': Solid concentration ratio <B>\theta<sub>s</sub></B>');
    FOM.result('pg6').set('titleparamindicator', false);
    FOM.result('pg6').set('xlabel', 'Normalized x-coordinate (unitless)');
    FOM.result('pg6').set('xlabelactive', true);
    FOM.result('pg6').set('ylabel', 'Normalized y-coordinate (unitless)');
    FOM.result('pg6').set('ylabelactive', true);
    FOM.result('pg6').set('showlegendsmaxmin', true);
    FOM.result('pg6').set('legendpos', 'rightdouble');
    FOM.result('pg6').feature('surf1').set('colortable', 'DiscoLight');
    FOM.result('pg6').feature('surf1').set('smooth', 'internal');
    FOM.result('pg6').feature('surf1').set('resolution', 'normal');

    FOM.result.create('pg7', 'PlotGroup1D');
    FOM.result('pg7').create('ptgr2', 'PointGraph');
    FOM.result('pg7').feature('ptgr2').set('data', 'dset1');
    FOM.result('pg7').feature('ptgr2').selection.set(4);
    FOM.result('pg7').feature('ptgr2').set('expr', 'phi_s - Iapp*RcFN(1,T)');    
    FOM.result('pg7').label('Cell voltage');
    FOM.result('pg7').set('titletype', 'manual');
    FOM.result('pg7').set('title', 'Cell voltage profile');
    FOM.result('pg7').set('xlabel', 'Time (s)');
    FOM.result('pg7').set('xlabelactive', true);
    FOM.result('pg7').set('ylabel', 'Cell voltage (V)');
    FOM.result('pg7').set('ylabelactive', true);
    FOM.result('pg7').feature('ptgr2').label('Plot of cell voltage');
    FOM.result('pg7').feature('ptgr2').set('descractive', true);
    FOM.result('pg7').feature('ptgr2').set('descr', 'Cell voltage');
    FOM.result('pg7').feature('ptgr2').set('titletype', 'custom');

    FOM.result.create('pg8', 'PlotGroup1D');
    FOM.result('pg8').create('lngr1', 'LineGraph');
    FOM.result('pg8').feature('lngr1').set('xdata', 'expr');
    FOM.result('pg8').feature('lngr1').selection.all;
    FOM.result('pg8').feature('lngr1').set('expr', 'thetass');
    FOM.result('pg8').label('theta_ss vs x');
    FOM.result('pg8').set('titletype', 'manual');
    FOM.result('pg8').set('title', 'Solid surface concentration ratio<B> \theta<sub>ss</sub></B>');
    FOM.result('pg8').set('xlabel', 'Normalized x-coordinate (unitless)');
    FOM.result('pg8').set('xlabelactive', true);
    FOM.result('pg8').set('ylabel', 'Concentration ratio (unitless)');
    FOM.result('pg8').set('ylabelactive', true);
    FOM.result('pg8').feature('lngr1').set('descr', 'thetass');
    FOM.result('pg8').feature('lngr1').set('titletype', 'manual');
    FOM.result('pg8').feature('lngr1').set('title', 'Solid surface concentration (mol/m^3)');
    FOM.result('pg8').feature('lngr1').set('xdataexpr', 'x/xnorm');
    FOM.result('pg8').feature('lngr1').set('xdataunit', '1');
    FOM.result('pg8').feature('lngr1').set('xdatadescractive', true);
    FOM.result('pg8').feature('lngr1').set('xdatadescr', 'x-coordinate (unitless)');
    FOM.result('pg8').feature('lngr1').set('smooth', 'none');
    FOM.result('pg8').feature('lngr1').set('resolution', 'normal');

    FOM.result.create('pg9', 'PlotGroup1D');
    FOM.result('pg9').create('lngr1', 'LineGraph');
    FOM.result('pg9').feature('lngr1').set('xdata', 'expr');
    FOM.result('pg9').feature('lngr1').selection.all;
    FOM.result('pg9').feature('lngr1').set('expr', 'if');
    FOM.result('pg9').label('if vs x');
    FOM.result('pg9').set('titletype', 'manual');
    FOM.result('pg9').set('title', 'Faradaic (intercalation) interface lithium flux <B>i<sub>f</sub></B>');
    FOM.result('pg9').set('xlabel', 'Normalized x-coordinate (unitless)');
    FOM.result('pg9').set('xlabelactive', true);
    FOM.result('pg9').set('ylabel', 'Flux (A)');
    FOM.result('pg9').set('ylabelactive', true);
    FOM.result('pg9').feature('lngr1').set('xdataexpr', 'x');
    FOM.result('pg9').feature('lngr1').set('xdataunit', 'm');
    FOM.result('pg9').feature('lngr1').set('xdatadescr', 'x-coordinate');
    FOM.result('pg9').feature('lngr1').set('smooth', 'none');
    FOM.result('pg9').feature('lngr1').set('resolution', 'normal');

    FOM.result.create('pg10', 'PlotGroup1D');
    FOM.result('pg10').create('lngr1', 'LineGraph');
    FOM.result('pg10').feature('lngr1').set('xdata', 'expr');
    FOM.result('pg10').feature('lngr1').selection.all;
    FOM.result('pg10').feature('lngr1').set('expr', 'idl');
    FOM.result('pg10').label('idl vs x');
    FOM.result('pg10').set('titletype', 'manual');
    FOM.result('pg10').set('title', 'Non-faradaic (double-layer) interface lithium flux <B>i<sub>dl</sub> </B>');
    FOM.result('pg10').set('xlabel', 'Normalized x-coordinate (unitless)');
    FOM.result('pg10').set('xlabelactive', true);
    FOM.result('pg10').set('ylabel', 'Flux (A)');
    FOM.result('pg10').set('ylabelactive', true);
    FOM.result('pg10').feature('lngr1').set('xdataexpr', 'x');
    FOM.result('pg10').feature('lngr1').set('xdataunit', 'm');
    FOM.result('pg10').feature('lngr1').set('xdatadescr', 'x-coordinate');
    FOM.result('pg10').feature('lngr1').set('smooth', 'none');
    FOM.result('pg10').feature('lngr1').set('resolution', 'normal');

    FOM.result.create('pg12', 'PlotGroup1D');
    FOM.result('pg12').create('ptgr2', 'PointGraph');
    FOM.result('pg12').feature('ptgr2').set('data', 'dset1');
    FOM.result('pg12').feature('ptgr2').selection.set(4);
    FOM.result('pg12').feature('ptgr2').set('expr', 'Iapp'); % 20200715
    FOM.result('pg12').label('Input current');
    FOM.result('pg12').set('titletype', 'manual');
    FOM.result('pg12').set('title', 'Input current profile');
    FOM.result('pg12').set('xlabel', 'Time (s)');
    FOM.result('pg12').set('xlabelactive', true);
    FOM.result('pg12').set('ylabel', 'Input current (A)');
    FOM.result('pg12').set('ylabelactive', true);
    FOM.result('pg12').feature('ptgr2').label('Plot of input current');
    FOM.result('pg12').feature('ptgr2').set('unit', '');
    FOM.result('pg12').feature('ptgr2').set('descractive', true);
    FOM.result('pg12').feature('ptgr2').set('descr', 'Input current');
    FOM.result('pg12').feature('ptgr2').set('titletype', 'custom');
    
    FOM.result.create('pg13', 'PlotGroup1D');
    FOM.result('pg13').create('lngr1', 'LineGraph');
    FOM.result('pg13').feature('lngr1').set('data', 'dset1');
    FOM.result('pg13').feature('lngr1').set('xdata', 'expr');
    FOM.result('pg13').feature('lngr1').selection.all;
    FOM.result('pg13').feature('lngr1').set('expr', 'phi_s-phi_e');
    FOM.result('pg13').label('phi_se vs x');
    FOM.result('pg13').set('titletype', 'manual');
    FOM.result('pg13').set('title', 'Phase potential difference <B> \phi<sub>s-e</sub></B>');
    FOM.result('pg13').set('xlabel', 'Normalized x-coordinate (unitless)');
    FOM.result('pg13').set('xlabelactive', true);
    FOM.result('pg13').set('ylabel', 'Potential (V)');
    FOM.result('pg13').set('ylabelactive', true);
    FOM.result('pg13').feature('lngr1').set('xdataexpr', 'x');
    FOM.result('pg13').feature('lngr1').set('xdataunit', 'm');
    FOM.result('pg13').feature('lngr1').set('xdatadescr', 'x-coordinate');
    FOM.result('pg13').feature('lngr1').set('smooth', 'none');
    FOM.result('pg13').feature('lngr1').set('resolution', 'normal');

    FOM.result.create('pg14', 'PlotGroup1D');
    FOM.result('pg14').create('lngr1', 'LineGraph');
    FOM.result('pg14').feature('lngr1').set('xdata', 'expr');
    FOM.result('pg14').feature('lngr1').selection.all;
    FOM.result('pg14').feature('lngr1').set('expr', 'thetasavg');
    FOM.result('pg14').label('theta_savg vs (x,t)');
    FOM.result('pg14').set('titletype', 'manual');
    FOM.result('pg14').set('title', 'Solid average concentration ratio<B> \theta<sub>s,avg</sub></B>');
    FOM.result('pg14').set('xlabel', 'Normalized x-coordinate (unitless)');
    FOM.result('pg14').set('xlabelactive', true);
    FOM.result('pg14').set('ylabel', 'Concentration ratio (unitless)');
    FOM.result('pg14').set('ylabelactive', true);
    FOM.result('pg14').feature('lngr1').set('descr', 'Average solid particle concentration ratio');
    FOM.result('pg14').feature('lngr1').set('xdataexpr', 'x');
    FOM.result('pg14').feature('lngr1').set('xdataunit', 'm');
    FOM.result('pg14').feature('lngr1').set('xdatadescr', 'x-coordinate');
    FOM.result('pg14').feature('lngr1').set('smooth', 'none');
    FOM.result('pg14').feature('lngr1').set('resolution', 'normal');

    FOM.result.create('pg15', 'PlotGroup2D');
    FOM.result('pg15').create('surf1', 'Surface');
    FOM.result('pg15').feature('surf1').set('expr', 'socavg');
    FOM.result('pg15').label('theta_savg vs (x,r,t)');
    FOM.result('pg15').set('titletype', 'custom');
    FOM.result('pg15').set('typeintitle', false);
    FOM.result('pg15').set('descriptionintitle', false);
    FOM.result('pg15').set('suffixintitle', ': Solid average concentration ratio <B>\theta<sub>s,avg</sub></B>');
    FOM.result('pg15').set('titleparamindicator', false);
    FOM.result('pg15').set('xlabel', 'Normalized x-coordinate (unitless)');
    FOM.result('pg15').set('xlabelactive', true);
    FOM.result('pg15').set('ylabel', 'Normalized y-coordinate (unitless)');
    FOM.result('pg15').set('ylabelactive', true);
    FOM.result('pg15').set('showlegendsmaxmin', true);
    FOM.result('pg15').set('legendpos', 'rightdouble');
    FOM.result('pg15').feature('surf1').set('colortable', 'DiscoLight');
    FOM.result('pg15').feature('surf1').set('smooth', 'internal');
    FOM.result('pg15').feature('surf1').set('resolution', 'normal');



    FOM.result.create('pg17', 'PlotGroup1D');    
    FOM.result('pg17').create('ptgr2', 'PointGraph');
    FOM.result('pg17').feature('ptgr2').set('data', 'dset1');
    FOM.result('pg17').feature('ptgr2').selection.set(3);
    FOM.result('pg17').feature('ptgr2').set('expr', 'thetasavg_pos');
    FOM.result('pg17').label('theta_savg_pos vs (t)');
    FOM.result('pg17').set('titletype', 'manual');
    FOM.result('pg17').set('title', 'Positive-electrode solid average concentration ratio<B> \theta<sub>s,avg</sub></B>');
    FOM.result('pg17').set('xlabel', 'Time (s)');
    FOM.result('pg17').set('xlabelactive', true);
    FOM.result('pg17').set('ylabel', 'Concentration ratio (unitless)');
    FOM.result('pg17').set('ylabelactive', true);
    FOM.result('pg17').feature('ptgr2').label('Plot of cathode average concentration ratio');
    FOM.result('pg17').feature('ptgr2').set('descr', 'thetasavg_pos');
    FOM.result('pg17').feature('ptgr2').set('titletype', 'custom');

    FOM.result.export.create('data2', 'Data');
    FOM.result.export('data2').label('phis');
    FOM.result.export('data2').set('expr', {'phi_s'});
    FOM.result.export('data2').set('unit', {'V'});
    FOM.result.export('data2').set('descr', {'Dependent variable phi_s'});
    FOM.result.export('data2').set('filename', sprintf('%s/Phis.txt',pwd));
    
    FOM.result.export.create('data3', 'Data');
    FOM.result.export('data3').label('phie');
    FOM.result.export('data3').set('expr', {'phi_e'});
    FOM.result.export('data3').set('unit', {'V'});
    FOM.result.export('data3').set('descr', {'Dependent variable phi_e'});
    FOM.result.export('data3').set('filename', sprintf('%s/Phie.txt',pwd));
    
    FOM.result.export.create('data4', 'Data');
    FOM.result.export('data4').label('thetass');
    FOM.result.export('data4').set('expr', {'thetass'});
    FOM.result.export('data4').set('unit', {'1'});
    FOM.result.export('data4').set('descr', {''});
    FOM.result.export('data4').set('filename', sprintf('%s/Thetass.txt',pwd));
    
    FOM.result.export.create('data5', 'Data');
    FOM.result.export('data5').label('thetae');
    FOM.result.export('data5').set('expr', {'theta_e'});
    FOM.result.export('data5').set('unit', {'1'});
    FOM.result.export('data5').set('descr', {''});
    FOM.result.export('data5').set('filename', sprintf('%s/Thetae.txt',pwd));
    
    FOM.result.export.create('data6', 'Data');
    FOM.result.export('data6').label('if');
    FOM.result.export('data6').set('expr', {'if'});
    FOM.result.export('data6').set('unit', {'A'});
    FOM.result.export('data6').set('descr', {'Dependent variable if'});
    FOM.result.export('data6').set('filename', sprintf('%s/If.txt',pwd));
    
    FOM.result.export.create('data7', 'Data');
    FOM.result.export('data7').label('phise');
    FOM.result.export('data7').set('expr', {'phi_s-phi_e'});
    FOM.result.export('data7').set('unit', {'V'});
    FOM.result.export('data7').set('descr', {''});
    FOM.result.export('data7').set('filename', sprintf('%s/Phise.txt',pwd));
    
    FOM.result.export.create('data8', 'Data');
    FOM.result.export('data8').label('idl');
    FOM.result.export('data8').set('expr', {'idl'});
    FOM.result.export('data8').set('unit', {'A'});
    FOM.result.export('data8').set('descr', {'Dependent variable idl'});
    FOM.result.export('data8').set('filename', sprintf('%s/Idl.txt',pwd));
    
    FOM.result.export.create('data9', 'Data');
    FOM.result.export('data9').label('Input current');
    FOM.result.export('data9').set('expr', {'Iapp'});
    FOM.result.export('data9').set('unit', {''});
    FOM.result.export('data9').set('descr', {'Constant current'});
    FOM.result.export('data9').set('filename', sprintf('%s/Iapp.txt',pwd));

    FOM.result.export.create('data10', 'Data');
    FOM.result.export('data10').label('thetasavg_neg vs (t)');
    FOM.result.export('data10').set('expr', {'thetasavg_neg'});
    FOM.result.export('data10').set('unit', {'1'});
    FOM.result.export('data10').set('descr', {''});
    FOM.result.export('data10').set('filename', sprintf('%s/ThetasAvgNeg.txt',pwd));
    
    FOM.result.export.create('data11', 'Data');
    FOM.result.export('data11').label('thetasavg_pos vs (t)');
    FOM.result.export('data11').set('expr', {'thetasavg_pos'});
    FOM.result.export('data11').set('unit', {'1'});
    FOM.result.export('data11').set('descr', {''});
    FOM.result.export('data11').set('filename', sprintf('%s/ThetasAvgPos.txt',pwd));
    
    FOM.result.export.create('data12', 'Data');
    FOM.result.export('data12').label('thetasavg');
    FOM.result.export('data12').set('expr', {'thetasavg'});
    FOM.result.export('data12').set('unit', {'1'});
    FOM.result.export('data12').set('descr', {''});
    FOM.result.export('data12').set('filename', sprintf('%s/ThetasAvg.txt',pwd));
  end

  function organize
    msg('\nOrganizing GUI...');
    FOM.nodeGroup.create('grp3', 'GlobalDefinitions');
    FOM.nodeGroup('grp3').set('type', 'func');
    FOM.nodeGroup('grp3').placeAfter([]);
    FOM.nodeGroup('grp3').label('Input-current profiles');
    FOM.nodeGroup('grp3').add('func', 'Ides');
    
    FOM.nodeGroup.create('grp3b', 'GlobalDefinitions');
    FOM.nodeGroup('grp3b').set('type', 'func');
    FOM.nodeGroup('grp3b').placeAfter([]);
    FOM.nodeGroup('grp3b').label('Temperature profiles');
    FOM.nodeGroup('grp3b').add('func', 'Temperature');

    FOM.nodeGroup.create('grp4', 'GlobalDefinitions');
    FOM.nodeGroup('grp4').set('type', 'func');
    FOM.nodeGroup.move('grp4', 1);
    FOM.nodeGroup('grp4').placeAfter([]);
    FOM.nodeGroup('grp4').label('Functions');
    FOM.nodeGroup('grp4').add('func', 'an1');
    FOM.nodeGroup('grp4').add('func', 'step1');
    
    FOM.nodeGroup.create('grp9', 'Definitions', 'mod1d');
    FOM.nodeGroup('grp9').set('type', 'cpl');
    FOM.nodeGroup('grp9').placeAfter('cpl', 'linext2');
    FOM.nodeGroup('grp9').label('Integration along x');
    FOM.nodeGroup('grp9').add('cpl', 'negint');
    FOM.nodeGroup('grp9').add('cpl', 'posint');
    
    FOM.nodeGroup.create('grp5', 'Definitions', 'mod2d');
    FOM.nodeGroup('grp5').set('type', 'cpl');
    FOM.nodeGroup('grp5').placeAfter('cpl', 'linext4');
    FOM.nodeGroup('grp5').label('Integration along y');
    FOM.nodeGroup('grp5').add('cpl', 'linproj2');
  
    FOM.nodeGroup('grp5').add('cpl', 'genproj2');
    
    FOM.nodeGroup.create('grp7', 'Results');
    FOM.nodeGroup('grp7').set('type', 'export');
    FOM.nodeGroup.move('grp7', 6);
    FOM.nodeGroup('grp7').placeAfter([]);
    FOM.nodeGroup('grp7').label('Data');
    FOM.nodeGroup('grp7').add('export', 'data9');
    FOM.nodeGroup('grp7').add('export', 'data6');
    FOM.nodeGroup('grp7').add('export', 'data8');
    FOM.nodeGroup('grp7').add('export', 'data2');
    FOM.nodeGroup('grp7').add('export', 'data3');
    FOM.nodeGroup('grp7').add('export', 'data7');
    FOM.nodeGroup('grp7').add('export', 'data5');
    FOM.nodeGroup('grp7').add('export', 'data4');
    FOM.nodeGroup('grp7').add('export', 'data10');
    FOM.nodeGroup('grp7').add('export', 'data11');
    FOM.nodeGroup('grp7').add('export', 'data12');
    
    FOM.nodeGroup.create('grp8', 'Results');
    FOM.nodeGroup('grp8').set('type', 'plotgroup');
    FOM.nodeGroup('grp8').placeAfter('plotgroup', 'pg8');
    FOM.nodeGroup('grp8').label('theta_savg');
    FOM.nodeGroup('grp8').add('plotgroup', 'pg14');
%     FOM.nodeGroup('grp8').add('plotgroup', 'pg16');
    FOM.nodeGroup('grp8').add('plotgroup', 'pg17');
    FOM.nodeGroup('grp8').add('plotgroup', 'pg15');  
    FOM.nodeGroup.move('grp3b', 2);
  end
  function msg(theText)
    if debugFlag
      fprintf(theText);
    end
  end
end

function [model, paramnamesflat] = explodeMSMR(model)
    %EXPLODEMSMR We need to handle MSMR models differently in
    % COMSOL. COMSOL does not support vector parameters, so this function
    % converts each length-J vector MSMR parameter to J scalar parameters.

    regs = {'pos'};
    params = {'U0','X','omega','k0','alpha'};
    paramnamesflat = struct;

    for r = 1:length(regs)
        reg = regs{r};

        % Build a set of scalar functions for each MSMR parameter.
        for k = 1:length(params)
            name = params{k};

            if model.MSMR
                % Extract the value and activation energy of the MSMR parameter 
                % into separate vectors.
                fcn = model.function.pos.(name);
                obj = reverseFcn(fcn);
                if ~strcmpi(obj.type,'vector') && ~strcmpi(obj.type,'scalar')
                    error('Failed to reverse MSMR parameter function %s.',name);
                end
                val = obj.value(:);
                Ea = obj.Ea(:);
                J = length(val);
                if length(Ea)==1
                    Ea = ones(J,1)*Ea;
                end
    
                % Construct a set of J scalar functions giving the value of each
                % element in the vector.
                paramnamesflat.(reg).(name) = cell(1,J);
                for i = 1:J
                    name_i = sprintf('%s%d',name,i);
                    fcn_i = str2func( ...
                        sprintf( ...
                            '@(x,T)(%g*exp(%g*(1/298.15-1/T)/8.31446))', ...
                            val(i),Ea(i) ...
                        ) ...
                    );
                    model.function.(reg).(name_i) = fcn_i;
                    paramnamesflat.(reg).(name){i} = name_i;
                end

                % Remove vector parameter (COMSOL does not support vectors!)
                model.function.(reg) = rmfield(model.function.pos,name);
            else
                paramnamesflat.(reg).(name) = {'null__'};
            end
        end

        % Store the number of MSMR galleries for convienent access.
        if model.MSMR
            model.function.(reg).J = str2func(sprintf('@(x,T)(%g)',J));
        end
    end
end