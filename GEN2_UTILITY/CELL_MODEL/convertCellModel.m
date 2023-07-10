function outModel = convertCellModel(inModel,outModelType)
%CONVERTCELLMODEL Transform cell model into another format.
%
% outModel = convertCellModel(inModel,outModelType) converts the cell
%   model INMODEL into the format specified by the OUTMODELTYPE string.
%
% Supported model formats:
%
%   'P2DM'    Psuedo-two-dimensional model
%   'WRM'     Warburg resistance model
%   'SWRM'    Spline WRM
%   'RLWRM'   Reduced-layer WRM
%   'RLSWRM'  Reduced-layer SWRM
%   'LLPM'    Legacy lumped-parameter model
%
% Note: this function can only down-convert models in this list;
% models cannot be up-converted due to irreversible parameter lumping.
%
% -- Changelog --
% 2023.06.22 | Created | Wesley Hileman <whileman@uccs.edu>

isModelType = @(x)any(strcmp(x,{'P2DM','WRM','SWRM','RLWRM','RLSWRM','LLPM'}));
parser = inputParser;
parser.addRequired('inStruct',@isstruct);
parser.addRequired('outModelType',isModelType);
parser.addParameter('MakeMSMROCP',true,@(x)isscalar(x)&&islogical(x));
parser.parse(inModel,outModelType);
arg = parser.Results; % structure of validated arguments

% Determine the type of model that was input.
inModelType = getCellModelType(inModel);

% No need to convert if model already in the requested form.
if strcmp(inModelType,outModelType)
    outModel = inModel;
    return;
end

% Map models to "reduction" factor; we can't perform
% the conversion when redfac.(outModelType) <= redfac.(inModelType).
redfac.P2DM = 0;
redfac.WRM = 1;
redfac.SWRM = 2;
redfac.RLWRM = 2;
redfac.RLSWRM = 3;
redfac.LLPM = 4;
redfacIN = redfac.(inModelType);
redfacOUT = redfac.(outModelType);
if redfacOUT <= redfacIN
    error('Upconversion: Cannot convert model type from %s to %s.', ...
        inModelType,outModelType);
end

% OK to convert.
if strcmp(inModelType,'P2DM') && strcmp(outModelType,'WRM')
    outModel = genWRM(inModel);
elseif strcmp(inModelType,'P2DM') && strcmp(outModelType,'RLWRM')
    tmpWRM = genWRM(inModel);
    outModel = genRLWRM(tmpWRM);
elseif strcmp(inModelType,'P2DM') && strcmp(outModelType,'SWRM')
    tmpWRM = genWRM(inModel);
    outModel = genSWRM(tmpWRM);
elseif strcmp(inModelType,'P2DM') && strcmp(outModelType,'RLSWRM')
    tmpWRM = genWRM(inModel);
    tmpRLWRM = genRLWRM(tmpWRM);
    outModel = genSWRM(tmpRLWRM);
elseif strcmp(inModelType,'WRM') && strcmp(outModelType,'RLWRM')
    outModel = genRLWRM(inModel);
elseif strcmp(inModelType,'WRM') && strcmp(outModelType,'SWRM')
    outModel = genSWRM(inModel);
elseif strcmp(inModelType,'WRM') && strcmp(outModelType,'RLSWRM')
    tmpRLWRM = genRLWRM(inModel);
    outModel = genSWRM(tmpRLWRM);
elseif strcmp(inModelType,'P2DM') && strcmp(outModelType,'LLPM')
    tmp = genWRM(inModel);
    outModel = genLLPM(tmp,arg);
elseif any(strcmp(inModelType,{'WRM','RLWRM'})) && strcmp(outModelType,'LLPM')
    outModel = genLLPM(inModel,arg);
else
    error('Not implemented: cannot convert %s to %s.', ...
        inModelType,outModelType);
end

outModel.type__ = 'cellModel';
outModel.origin__ = 'convertCellModel';
outModel.arg__ = arg;

end

function WRM = genWRM(P2DM)
%GENWORM Convert P2DM to WRM.

R = TB.const.R;
F = TB.const.F;

WRM = struct;
WRM.metadata = P2DM.metadata;
WRM.metadata.cell.type = 'WRM';
WRM.metadata.cell.lumpedParams = true;
WRM.parameters = struct;

p = getNumericParams(P2DM.parameters);
secNames = fieldnames(P2DM.metadata.section);
for s = 1:length(secNames)
    secName = secNames{s};
    secMeta = P2DM.metadata.section.(secName);

    l = struct;
    switch secMeta.type
        case 'Package'
            l.R0 = genNumericParam('R0',p.(secName).R0,0,'Ohm');
            l.L0 = genNumericParam('L0',p.(secName).L0,0,'H');
        case 'Global'
            kD = 2*R*(p.const.t0plus-1)*(1+p.const.dlnfdlnc)/F; % [V/K]
            psi = F*p.const.De*p.const.ce0/p.const.kappa/(1-p.const.t0plus)/p.const.Tref; % [V/K]
            l.Q = genNumericParam('Q',p.pos.sEps*p.const.A*p.pos.L*F*p.pos.csmax*abs(p.pos.theta100-p.pos.theta0)/3600,0,'A*h');
            l.psi = genNumericParam('psi',psi,0,'V/K');
            l.W = genNumericParam('W',-kD/psi,0,'unitless');
            l.Tref = genNumericParam('Tref',p.const.Tref,0,'K');
        case 'Electrode3D'
            as = 3*p.(secName).sEps/p.(secName).Rs;
            l.sigma = genNumericParam('sigma',p.const.A*p.(secName).sigma*(p.(secName).sEps)^p.(secName).brugSigma/p.(secName).L,0,'1/Ohm');
            l.kappa = genNumericParam('kappa',p.const.A*p.const.kappa*(p.(secName).eEps)^p.(secName).brugDeKappa/p.(secName).L,0,'1/Ohm');
            l.tauW = genNumericParam('tauW',p.(secName).eEps*(p.(secName).L)^2/p.const.De/(p.(secName).eEps)^p.(secName).brugDeKappa,0,'s');
            l.Dsref = genNumericParam('Dsref',p.(secName).Dsref/p.(secName).Rs^2,0,'1/s');
            l.nF = genNumericParam('nF',p.(secName).nF,0,'unitless');
            l.nDL = genNumericParam('nDL',p.(secName).nDL,0,'unitless');
            l.Rf = genNumericParam('Rf',p.(secName).Rf/as/p.const.A/p.(secName).L,0,'Ohm');
            l.Rdl = genNumericParam('Rdl',p.(secName).Rdl/as/p.const.A/p.(secName).L,0,'Ohm');
            l.Cdl = genNumericParam('Cdl',p.(secName).Cdl*as*p.const.A*p.(secName).L,0,'F');
            l.alpha = genNumericParam('alpha',p.(secName).alpha,0,'unitless');
            l.theta0 = genNumericParam('theta0',p.(secName).theta0,0,'unitless');
            l.theta100 = genNumericParam('theta100',p.(secName).theta100,0,'unitless');
            l.k0 = genNumericParam('k0',p.const.A*p.(secName).L*as*p.(secName).k0,0,'A');
        case 'Electrode2D'
            l.k0 = genNumericParam('k0',p.(secName).gamma*p.const.A*F*p.(secName).k0*(p.(secName).cs0)^p.(secName).alpha*p.const.ce0^(1-p.(secName).alpha),0,'A');
            l.alpha = genNumericParam('alpha',p.(secName).alpha,0,'unitless');
            l.nDL = genNumericParam('nDL',p.(secName).nDL,0,'unitless');
            l.Rf = genNumericParam('Rf',p.(secName).Rf/p.(secName).gamma/p.const.A,0,'Ohm');
            l.Rdl = genNumericParam('Rdl',p.(secName).Rdl/p.(secName).gamma/p.const.A,0,'Ohm');
            l.Cdl = genNumericParam('Cdl',p.(secName).Cdl*p.(secName).gamma*p.const.A,0,'F');
        case 'ElectrolyteLayer'
            l.kappa = genNumericParam('kappa',p.const.A*p.const.kappa*(p.(secName).eEps)^p.(secName).brugDeKappa/p.(secName).L,0,'1/Ohm');
            l.tauW = genNumericParam('tauW',p.(secName).eEps*(p.(secName).L)^2/p.const.De/(p.(secName).eEps)^p.(secName).brugDeKappa,0,'s');
    end % switch
    WRM.parameters.(secName) = l;

    if strcmp(secMeta.type,'Electrode3D')
        % Copy OCP.
        if isfield(P2DM.parameters.(secName),'Uocv')
            WRM.parameters.(secName).Uocv = P2DM.parameters.(secName).Uocv;
        end
        if isfield(P2DM.parameters.(secName),'dUocv')
            WRM.parameters.(secName).dUocv = P2DM.parameters.(secName).dUocv;
        end
        if isfield(P2DM.parameters.(secName),'d2Uocv')
            WRM.parameters.(secName).d2Uocv = P2DM.parameters.(secName).d2Uocv;
        end
        if isfield(P2DM.parameters.(secName),'U0')
            WRM.parameters.(secName).U0 = P2DM.parameters.(secName).U0;
        end
        if isfield(P2DM.parameters.(secName),'X')
            WRM.parameters.(secName).X = P2DM.parameters.(secName).X;
        end
        if isfield(P2DM.parameters.(secName),'omega')
            WRM.parameters.(secName).omega = P2DM.parameters.(secName).omega;
        end
    end
end % for sec

end

function RLWRM = genRLWRM(WRM)
%GENRLWORM Convert WRM to RLWRM (Reduced-Layer WRM).

RLWRM = WRM;
RLWRM.metadata.cell.type = 'RLWRM';

% Combine dll/sep layers into single eff layer.
p = getNumericParams(RLWRM.parameters);
k1 = 1+p.dll.kappa/p.sep.kappa;
k2 = 1+p.sep.kappa/p.dll.kappa;
eff.tauW = genNumericParam('tauW',p.dll.tauW*k1+p.sep.tauW*k2,0);
eff.kappa = genNumericParam('kappa', ...
    p.dll.kappa*p.sep.kappa/(p.dll.kappa+p.sep.kappa),0);
RLWRM.parameters.sep = eff;
RLWRM.parameters = renameStructField(RLWRM.parameters,'sep','eff');
RLWRM.parameters = rmfield(RLWRM.parameters,'dll');

% Update metadata.
RLWRM.metadata.section = renameStructField( ...
    RLWRM.metadata.section,'sep','eff');
RLWRM.metadata.section = rmfield(RLWRM.metadata.section,'dll');
RLWRM.metadata.section.eff.name = 'eff';

end

function SWRM = genSWRM(WRM)
%GENSWRM Convert WRM (or RLWRM) to SWRM (or RLSWRM).

SWRM = WRM;
if strcmpi(WRM.metadata.cell.type,'WRM')
    SWRM.metadata.cell.type = 'SWRM';
elseif strcmpi(WRM.metadata.cell.type,'RLWRM')
    SWRM.metadata.cell.type = 'RLSWRM';
else
    error('Invalid cell model type. Must be WRM or RLWRM.');
end

params = getNumericParams(WRM.parameters);
if isfield(params.const,'Tref')
    TrefdegC = params.const.Tref-273.15;
else
    TrefdegC = 25;
end
secNames = fieldnames(WRM.metadata.section);
for s = 1:length(secNames)
    secName = secNames{s};
    secMeta = WRM.metadata.section.(secName);

    if strcmp(secMeta.type,'Electrode3D')
        if ~strcmpi(secMeta.ocp.type,'msmr') || ...
           ~strcmpi(secMeta.kinetics,'msmr') || ...
           ~strcmpi(secMeta.solidDiffusion,'msmr')
            error(['Not implemented: cannot generate spline model from ' ...
                'non-MSMR model (OCP, kinetics, and Ds must be MSMR).'])
        end

        electrode = MSMR(params.(secName));
        thetaSpline = linspace( ...
            electrode.zmin,electrode.zmax,electrode.J).';

        % Convert MSMR kinetics.
        % Add spline kinetics parameters.
        ctData = electrode.Rct(params.(secName), ...
            'theta',thetaSpline,'TdegC',TrefdegC);
        k0Spline = ctData.i0;
        SWRM.parameters.(secName).k0SplineTheta = genNumericParam( ...
            'k0SplineTheta',thetaSpline(:),0,'unitless');
        SWRM.parameters.(secName).k0Spline = genNumericParam( ...
            'k0Spline',k0Spline(:),0,'A');
        % Remove MSMR kinetics parameters.
        SWRM.parameters.(secName) = rmfield( ...
            SWRM.parameters.(secName),'k0');

        % Convert MSMR solid diffusivity.
        % Add spline diffusivity parameters.
        dsData = electrode.Ds(params.(secName), ...
            'theta',thetaSpline(:),'TdegC',TrefdegC);
        DsSpline = dsData.Ds;
        SWRM.parameters.(secName).DsSplineTheta = genNumericParam( ...
            'DsSplineTheta',thetaSpline(:),0,'unitless');
        SWRM.parameters.(secName).DsSpline = genNumericParam( ...
            'DsSpline',DsSpline(:),0,'1/s');
        % Remove MSMR diffusivity parameters.
        SWRM.parameters.(secName) = rmfield( ...
            SWRM.parameters.(secName),'Dsref');

        % Update metadata.
        SWRM.metadata.section.(secName).kinetics = 'spline';
        SWRM.metadata.section.(secName).solidDiffusion = 'spline';
    end
end % for

end

function paramStruct = genNumericParam(name,value,Eact,unit)
%GENNUMERICPARAM Create parameter structure for numeric data.

if ~exist('unit','var')
    unit = '';
end

paramStruct.Name = name;
paramStruct.Value.type = 'Numeric';
paramStruct.Value.isScalar = length(value)==1;
paramStruct.Value.value = value;
paramStruct.Eact = Eact;
paramStruct.Unit = unit;
if isfield(TB.const.LUMPED_PARAMETER_METADATA,name)
    paramMeta = TB.const.LUMPED_PARAMETER_METADATA.(name);
    paramStruct.Unit = paramMeta.u; % we know better than the user
    paramStruct.Description = paramMeta.d;
    paramStruct.Latex = paramMeta.l;
end

end

function numStruct = getNumericParams(paramsStruct)
%GETNUMERICPARAMS Extract values of numeric parameters from cell model.

secNames = fieldnames(paramsStruct);
for s = 1:length(secNames)
    secName = secNames{s};
    secStruct = paramsStruct.(secName);
    paramNames = fieldnames(secStruct);
    for p = 1:length(paramNames)
        paramName = paramNames{p};
        paramStruct = secStruct.(paramName);
        if strcmp(paramStruct.Value.type,'Numeric')
            numStruct.(secName).(paramName) = paramStruct.Value.value;
        end
    end % for
end % for

end

function LLPM = genLLPM(WORM,arg)
%GENLLPM Convert WORM or RLWORM to LLPM.

cellModel = WORM;
if strcmp(WORM.metadata.cell.type,'RLWRM')
    % Expand eff into dll and sep for legacy implemtation
    % of genFOM. The layers will have identical parameters
    % to emulate a single eff layer.
    p = getNumericParams(cellModel.parameters);
    tauW = genNumericParam('tauW',p.eff.tauW/4,0);
    kappa = genNumericParam('kappa',p.eff.kappa*2,0);
    cellModel.parameters = rmfield(cellModel.parameters,'eff');
    cellModel.parameters.dll = struct('tauW',tauW,'kappa',kappa);
    cellModel.parameters.sep = struct('tauW',tauW,'kappa',kappa);
    cellModel.parameters = orderfields( ...
        cellModel.parameters,{'const','neg','dll','sep','pos'});
    % Update metadata.
    cellModel.metadata.section.dll = cellModel.metadata.section.eff;
    cellModel.metadata.section.sep = cellModel.metadata.section.eff;
    cellModel.metadata.section = rmfield(cellModel.metadata.section,'eff');
    cellModel.metadata.section = orderfields( ...
        cellModel.metadata.section,{'const','neg','dll','sep','pos'});
end % if

LLPM.const.R = TB.const.R;
LLPM.const.F = TB.const.F;
LLPM.const.T = cellModel.parameters.const.Tref.Value.value;
LLPM.name = cellModel.metadata.cell.name;
LLPM.lumped = cellModel.metadata.cell.lumpedParams;
LLPM.MSMR = strcmpi(cellModel.metadata.section.pos.ocp.type,'MSMR');
LLPM.TM = false;

% Convert each section.
secNames = fieldnames(cellModel.metadata.section);
for s = 1:length(secNames)
    secName = secNames{s};
    secMeta = cellModel.metadata.section.(secName);
    sec = cellModel.parameters.(secName);
    paramNames = fieldnames(sec);
    for p = 1:length(paramNames)
        paramName = paramNames{p};
        param = sec.(paramName);
        switch(param.Value.type)
            case 'Numeric'
                LLPM.function.(secName).(paramName) = ... 
                    numericParam2function(param,LLPM.const.T);
            case 'LUT'
                error(['Not implemented: %s.%s\n' ...
                    'LUT value not yet implmented.'], ...
                    secName,paramName);
            case 'Function'
                error(['Not implemented: %s.%s\n' ...
                    'Function value not yet implemented.'], ...
                    secName,paramName);
            otherwise
                warning(['Unrecognized value type: %s.%s.\n' ...
                    'Type string: %s\n' ...
                    'Skipping conversion.'], ...
                    secName,paramName,param.Value.type);
        end % switch
    end % for parameter

    % Additional processing for porous electrodes.
    if strcmp(secMeta.type,'Electrode3D')
        % Add SOC function to electrode
        theta0 = LLPM.function.(secName).theta0(); 
        theta100 = LLPM.function.(secName).theta100();
        LLPM.function.(secName).soc = str2func( ...
            sprintf('@(x,T)(%g+x*(%g))',...
            theta0,theta100-theta0));

        if LLPM.MSMR && arg.MakeMSMROCP
            % Compute OCP functions from MSMR parameters.
            [UocpFcn, dUocpFcn] = makeOCP(LLPM.function.(secName));
            LLPM.function.(secName).Uocp = UocpFcn;
            LLPM.function.(secName).dUocp = dUocpFcn;
        end % if
    end % if
end % for section

% Make cell-level OCP function.
if isfield(LLPM.function.pos,'Uocp')
    LLPM = makeOCV(LLPM);
end

end % genLLPM()

function fcn = numericParam2function(param,Tref)
%NUMERIC2FUNCTION Convert numeric parameter to function.

value = param.Value.value;
if isscalar(value)
    valueString = sprintf('(%g)',value);
else
    valueString = ['([' sprintf('%g;',value) '])'];
end
Eact = param.Eact;
if all(Eact==0)
    EactString = '';
elseif isscalar(Eact)
    EactString = sprintf('%g',Eact);
else
    EactString = ['[' sprintf('%g;',Eact) ']'];
end
if isempty(EactString)
    fcn = str2func(sprintf('@(x,T)%s', valueString));
else
    fcn = str2func( ...
        sprintf('@(x,T)%s.*exp(%s*(1/%g-1/T)/8.3144621)', ...
        valueString,EactString,Tref) ...
    );
end

end % numeric2function()

function [Uocp, dUocp] = makeOCP(sec)
%MAKEDOCP Make electrode-level OCP function from MSMR parameter values.

electrode = MSMR(sec);
SOCvec = 0.001:0.0001:0.999;  % (stoichiometry, theta)
ocpData = electrode.ocp('theta',SOCvec,'TdegC',25);
OCVvec1 = ocpData.Uocp;
dOCVvec = ocpData.dUocp;
ocpData = electrode.ocp('theta',SOCvec,'TdegC',26);
OCVvec2 = ocpData.Uocp;

SOCstr = ['[',sprintf('%g;',SOCvec),']'];
OCVstr1 = ['[',sprintf('%g;',OCVvec1),']'];
OCVstr2 = ['[',sprintf('%g;',OCVvec2-OCVvec1),']'];
dOCVstr = ['[',sprintf('%g;',dOCVvec),']'];

% Make Uocp function, adding in temperature-dependence.
UocpFn = sprintf( ...
    '@(x,T) interp1(%s,%s+(T(:)-298.15).*%s,x,''linear'',''extrap'')', ...
     SOCstr,OCVstr1,OCVstr2);
Uocp = str2func(UocpFn);

% Make dUocp function, adding in temperature-dependence.
dUocpFn = sprintf( ...
    '@(x,T) interp1(%s,T(:).*%s,x,''linear'',''extrap'')',SOCstr,dOCVstr);
dUocp = str2func(dUocpFn); 

end

% Make cell-level OCV function from electrode-level OCP functions
function cellParams = makeOCV(cellParams)
    SOC = 0:0.001:1;
    posTheta = cellParams.function.pos.soc(SOC,298.15);
    posOCP = cellParams.function.pos.Uocp(posTheta,298.15);
    cellOCV = posOCP ;
    SOCstr = sprintf(' %g;',SOC);
    SOCstr = SOCstr(1:end-1); % remove trailing ';'
    OCVstr = sprintf(' %g;',cellOCV);
    OCVstr = OCVstr(1:end-1); % ditto
    OCVfn = sprintf( ...
        '@(x,T) interp1([%s],[%s],x,''linear'',''extrap'')',SOCstr,OCVstr);
    cellParams.function.const.Uocv = str2func(OCVfn);
end