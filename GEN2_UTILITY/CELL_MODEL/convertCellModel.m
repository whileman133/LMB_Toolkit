function outModel = convertCellModel(inModel,outModelType,varargin)
%CONVERTCELLMODEL Transform cell model into another format.
%
% outModel = convertCellModel(inModel,outModelType) converts the cell
%   model INMODEL into the format specified by the OUTMODELTYPE string.
%
% outModel = convertCellModel(inModel,'SolidDiffusionModel',dStruct)
%   converts the solid diffusion model of inModel into the format specified
%   by DSTRUCT. DSTRUCT is a structure with the `type` field and any other
%   fields needed by the specific format of model requested. If all other
%   fields are optional, then DSTRUCT may be a character vector specifying
%   the model format.
%
% outModel = convertCellModel(inModel,'KineticsModel',kStruct)
%   converts the kinetics model of inModel into the format specified
%   by KSTRUCT. KSTRUCT is a structure with the `type` field and any other
%   fields needed by the specific format of model requested. If all other
%   fields are optional, then KSTRUCT may be a character vector specifying
%   the model format.
%
% Supported model formats:
%
%   'P2DM'    Psuedo-two-dimensional model
%   'WRM'     Warburg resistance model
%   'RLWRM'   Reduced-layer WRM
%   'LLPM'    Legacy lumped-parameter model
%
% Supported solid-diffusion and kinetics model formats:
%
%   'MSMR'    Baker-Verbrugge MSMR model.
%   'linear'  Log linear interpolation over lithiation.
%             Optional: lithiation vector `theta`
%   'spline'  Log cubic spline interpolation over lithiation.
%             Optional: lithiation vector `theta`
%   'gp'      Log Gaussian Process (GP) regression over lithiation.
%
% Note: this function can only down-convert models in this list;
% models cannot be up-converted due to irreversible parameter lumping.
%
% -- Changelog --
% 2023.09.28 | Scaling for lumped MSMR constants k0 | Wes H.
% 2023.07.24 | Add diffusion/kinetics model conversion | Wesley Hileman
% 2023.06.22 | Created | Wesley Hileman <whileman@uccs.edu>

isModelType = @(x)any(strcmp(x,{'P2DM','WRM','RLWRM','LLPM'}));
isDiffusionType = @(x)any(strcmpi(x,{'msmr','linear','spline','gp'}));
isKineticsType = @(x)any(strcmpi(x,{'msmr','linear','spline','gp'}));
parser = inputParser;
parser.addRequired('inStruct',@isstruct);
parser.addOptional('outModelType',[],isModelType);
parser.addParameter('SolidDiffusionModel',[],@(x)isstruct(x)||ischar(x));
parser.addParameter('KineticsModel',[],@(x)isstruct(x)||ischar(x));
parser.addParameter('MakeMSMROCP',true,@(x)isscalar(x)&&islogical(x));
parser.addParameter('LegacyExpandEff',true,@(x)isscalar(x)&&islogical(x));
parser.parse(inModel,outModelType,varargin{:});
arg = parser.Results; % structure of validated arguments

% Determine the type of model that was input.
inModelType = getCellModelType(inModel);

% Validate diffusion/kinetics parameters.
if strcmp(inModelType,'LLPM') && (...
    ~isempty(arg.SolidDiffusionModel) || ~isempty(arg.KineticsModel))
    error('Not implemented: cannot convert kinetics of LLPM.');
end
if ischar(arg.SolidDiffusionModel)
    arg.SolidDiffusionModel = struct('type',arg.SolidDiffusionModel);
end
if ischar(arg.KineticsModel)
    arg.KineticsModel = struct('type',arg.KineticsModel);
end
if ~isempty(arg.SolidDiffusionModel) && ~isDiffusionType(arg.SolidDiffusionModel.type)
    error('Invalid solid diffusion model');
end
if ~isempty(arg.KineticsModel) && ~isKineticsType(arg.KineticsModel.type)
    error('Invalid kinetics model.');
end

% 1. Convert the model.
if isempty(outModelType) || strcmp(inModelType,outModelType)
    % No need to convert if model already in the requested form.
    outModel = inModel;
else
    % Map models to "reduction" factor; we can't perform
    % the conversion when redfac.(outModelType) <= redfac.(inModelType).
    redfac.P2DM = 0;
    redfac.RLP2DM = 1;
    redfac.WRM = 2;
    redfac.RLWRM = 3;
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
    elseif strcmp(inModelType,'RLP2DM') && strcmp(outModelType,'RLWRM')
        outModel = genWRM(inModel);
    elseif strcmp(inModelType,'P2DM') && strcmp(outModelType,'RLWRM')
        tmpWRM = genWRM(inModel);
        outModel = genRLWRM(tmpWRM);
    elseif strcmp(inModelType,'WRM') && strcmp(outModelType,'RLWRM')
        outModel = genRLWRM(inModel);
    elseif any(strcmp(inModelType,{'P2DM','RLP2DM'})) && strcmp(outModelType,'LLPM')
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

% 2. Convert the diffusion/kinetics models.
outModel = convertSubModels(outModel,arg.SolidDiffusionModel,arg.KineticsModel);

end

function WRM = genWRM(P2DM)
%GENWORM Convert P2DM to WRM, or RLP2DM to RLWRM.

R = TB.const.R;
F = TB.const.F;

WRM = struct;
WRM.metadata = P2DM.metadata;
if isfield(P2DM.parameters,'eff')
    WRM.metadata.cell.type = 'RLWRM';
else
    WRM.metadata.cell.type = 'WRM';
end
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
            l.Vmin = genNumericParam('Vmin',p.const.Vmin,0,'V');
            l.Vmax = genNumericParam('Vmax',p.const.Vmax,0,'V');
            l.Tref = genNumericParam('Tref',p.const.Tref,0,'K');
        case 'Electrode3D'
            as = 3*p.(secName).sEps/p.(secName).Rs;
            l.sigma = genNumericParam('sigma',p.const.A*p.(secName).sigma*(p.(secName).sEps)^p.(secName).brugSigma/p.(secName).L,0,'1/Ohm');
            l.kappa = genNumericParam('kappa',p.const.A*p.const.kappa*(p.(secName).eEps)^p.(secName).brugDeKappa/p.(secName).L,0,'1/Ohm');
            l.tauW = genNumericParam('tauW',p.(secName).eEps*(p.(secName).L)^2/p.const.De/(p.(secName).eEps)^p.(secName).brugDeKappa,0,'s');
            l.Dsref = genNumericParam('Dsref',p.(secName).Dsref/p.(secName).Rs^2,0,'1/s');
            l.nF = genNumericParam('nF',p.(secName).nF,0,'unitless');
            l.tauF = genNumericParam('tauF',p.(secName).tauF,0,'1/s');
            l.nDL = genNumericParam('nDL',p.(secName).nDL,0,'unitless');
            l.tauDL = genNumericParam('tauDL',p.(secName).tauDL,0,'s');
            l.Rf = genNumericParam('Rf',p.(secName).Rf/as/p.const.A/p.(secName).L,0,'Ohm');
            l.Rdl = genNumericParam('Rdl',p.(secName).Rdl/as/p.const.A/p.(secName).L,0,'Ohm');
            l.Cdl = genNumericParam('Cdl',p.(secName).Cdl*as*p.const.A*p.(secName).L,0,'F');
            l.alpha = genNumericParam('alpha',p.(secName).alpha,0,'unitless');
            l.theta0 = genNumericParam('theta0',p.(secName).theta0,0,'unitless');
            l.theta100 = genNumericParam('theta100',p.(secName).theta100,0,'unitless');
            % Rate constant(s).
            k0 = p.const.A*p.(secName).L*as*p.(secName).k0;
            if all(isfield(p.(secName),{'X','omega'}))
                % MSMR model; normalize k0.
                k0 = k0.*(p.(secName).X/2).^(p.(secName).omega);
            end
            l.k0 = genNumericParam('k0',k0,0,'A');
        case 'Electrode2D'
            l.k0 = genNumericParam('k0',p.(secName).gamma*p.const.A*F*p.(secName).k0*(p.(secName).cs0)^p.(secName).alpha*p.const.ce0^(1-p.(secName).alpha),0,'A');
            l.alpha = genNumericParam('alpha',p.(secName).alpha,0,'unitless');
            l.nDL = genNumericParam('nDL',p.(secName).nDL,0,'unitless');
            l.tauDL = genNumericParam('tauDL',p.(secName).tauDL,0,'s');
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
        if isfield(P2DM.parameters.(secName),'Uocp')
            WRM.parameters.(secName).Uocp = P2DM.parameters.(secName).Uocp;
        end
        if isfield(P2DM.parameters.(secName),'dUocp')
            WRM.parameters.(secName).dUocp = P2DM.parameters.(secName).dUocp;
        end
        if isfield(P2DM.parameters.(secName),'d2Uocp')
            WRM.parameters.(secName).d2Uocp = P2DM.parameters.(secName).d2Uocp;
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

function outModel = convertSubModels(inModel,solidDiffusionModel,kineticsModel)
%CONVERTSUBMODELS Convert solid diffusion and/or kinetics models.

outModel = inModel;
if isempty(solidDiffusionModel) && isempty(kineticsModel)
    % No need to convert.
    return;
end

params = getNumericParams(inModel.parameters);
if isfield(params.const,'Tref')
    TrefdegC = params.const.Tref-273.15;
else
    TrefdegC = 25;
end
secNames = fieldnames(inModel.metadata.section);
for s = 1:length(secNames)
    secName = secNames{s};
    secMeta = inModel.metadata.section.(secName);

    if strcmp(secMeta.type,'Electrode3D')
        if ~strcmpi(secMeta.ocp.type,'msmr')
            error(['Not implemented: cannot generate diffusion/kinetics ' ...
                'model from non-MSMR OCP model (OCP must be MSMR).'])
        end

        electrode = MSMR(params.(secName));

        % Convert kinetics.
        if isempty(kineticsModel) || strcmpi(secMeta.kinetics.type,kineticsModel.type)
            % Skip.
        elseif any(strcmpi(kineticsModel.type,{'linear','spline'}))
            if isfield(kineticsModel,'theta') && ~isempty(kineticsModel.theta)
                theta = kineticsModel.theta(:);
            else
                theta = linspace( ...
                    electrode.zmin,electrode.zmax,electrode.J).';
            end
            ctData = electrode.Rct(params.(secName), ...
                'theta',theta(:),'TdegC',TrefdegC);
            k0 = ctData.i0;
            alpha = 0.5*ones(size(theta(:)));
            % Remove any old kinetics parameters.
            outModel.parameters.(secName) = rmfield( ...
                outModel.parameters.(secName), ...
                intersect( ...
                    {'k0','k0Linear','k0Spline','k0Theta','alphaLinear', ...
                    'alphaSpline'}, ...
                    fieldnames(outModel.parameters.(secName)) ...
                ) ...
            );
            % Add new kinetics parameters.
            outModel.parameters.(secName).k0Theta = genNumericParam( ...
                    'k0Theta',theta(:),0,'unitless');
            if strcmpi(kineticsModel.type,'linear')
                outModel.parameters.(secName).k0Linear = genNumericParam( ...
                    'k0Linear',k0(:),0,'A');
                outModel.parameters.(secName).alphaLinear = genNumericParam( ...
                    'alphaLinear',alpha(:),0,'unitless');
            else
                outModel.parameters.(secName).k0Spline = genNumericParam( ...
                    'k0Spline',k0(:),0,'A');
                outModel.parameters.(secName).alphaSpline = genNumericParam( ...
                    'alphaSpline',alpha(:),0,'unitless');
            end
            % Update metadata.
            outModel.metadata.section.(secName).kinetics = struct;
            outModel.metadata.section.(secName).kinetics.type = kineticsModel.type;
        else
            error('Not implemented: cannot generate %s kinetics model.', ...
                kineticsModel.type);
        end % if kinetics

        % Convert solid diffusion.
        if isempty(solidDiffusionModel) || strcmpi(secMeta.solidDiffusion.type,solidDiffusionModel.type)
            % Skip.
        elseif any(strcmpi(solidDiffusionModel.type,{'linear','spline'}))
            if isfield(solidDiffusionModel,'theta') && ~isempty(solidDiffusionModel.theta)
                theta = solidDiffusionModel.theta(:);
            else
                theta = linspace( ...
                    electrode.zmin,electrode.zmax,electrode.J).';
            end
            dsData = electrode.Ds(params.(secName), ...
                'theta',theta(:),'TdegC',TrefdegC);
            Ds = dsData.Ds;
            % Remove any old diffusivity parameters.
            outModel.parameters.(secName) = rmfield( ...
                outModel.parameters.(secName), ...
                intersect({'Dsref','DsLinear','DsSpline','DsTheta'},fieldnames(outModel.parameters.(secName))));
            % Add new diffusivity parameters.
            outModel.parameters.(secName).DsTheta = genNumericParam( ...
                    'DsTheta',theta(:),0,'unitless');
            if strcmpi(solidDiffusionModel.type,'linear')
                outModel.parameters.(secName).DsLinear = genNumericParam( ...
                    'DsLinear',Ds(:),0,'1/s');
            else
                outModel.parameters.(secName).DsSpline = genNumericParam( ...
                    'DsSpline',Ds(:),0,'1/s');
            end
            % Update metadata.
            outModel.metadata.section.(secName).solidDiffusion = struct;
            outModel.metadata.section.(secName).solidDiffusion.type = solidDiffusionModel.type;
        else
            error('Not implemented: cannot generate %s solid diffusion model.', ...
                solidDiffusionModel.type);
        end % if diffusion
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
if any(strcmp(WORM.metadata.cell.type,{'RLWRM','RLP2DM'})) && arg.LegacyExpandEff
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
        cellModel.parameters,{'const','neg','dll','sep','pos','pkg'});
    % Update metadata.
    cellModel.metadata.section.dll = cellModel.metadata.section.eff;
    cellModel.metadata.section.sep = cellModel.metadata.section.eff;
    cellModel.metadata.section = rmfield(cellModel.metadata.section,'eff');
    cellModel.metadata.section = orderfields( ...
        cellModel.metadata.section,{'const','neg','dll','sep','pos','pkg'});
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
                LLPM.function.(secName).(paramName) = ...
                    lut2function(param);
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

% Compute legacy parameters qe, kD from new parameters W, tauW.
T = LLPM.const.T;
W = LLPM.function.const.W(1,T);
psi = LLPM.function.const.psi(1,T);
kD = -psi*W;
LLPM.function.const.kD = str2func(sprintf('@(x,T)(%g)',kD));
for s = 1:length(secNames)
    secName = secNames{s};
    secMeta = cellModel.metadata.section.(secName);
    if any(strcmp(secMeta.type,{'Electrode3D','ElectrolyteLayer'}))
        kappa = LLPM.function.(secName).kappa(1,T);
        tau = LLPM.function.(secName).tauW(1,T);
        qe = psi*T*kappa*tau/3600;
        LLPM.function.(secName).qe = str2func(sprintf('@(x,T)(%g)',qe));
    end % for
end % for
% Legacy tab resistance parameter.
LLPM.function.const.Rc = LLPM.function.pkg.R0;

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

function fcn = lut2function(param)
%LUT2FUNCTION Convert lookup table to function.

str = sprintf( ...
    '@(x,T) (interp1(%s,%s,x,''linear'',''extrap''))',...
    mat2str(param.Value.x),mat2str(param.Value.y));
fcn = str2func(str);

end

function [Uocp, dUocp] = makeOCP(sec)
%MAKEDOCP Make electrode-level OCP function from MSMR parameter values.

electrode = MSMR(sec);
SOCvec = 0.001:0.0001:0.999;  % (stoichiometry, theta)
ocpData = electrode.ocp('theta',SOCvec,'TdegC',25);
OCVvec1 = ocpData.Uocp;
dOCVvec = ocpData.dUocp*ocpData.f/(ocpData.F/ocpData.R);
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