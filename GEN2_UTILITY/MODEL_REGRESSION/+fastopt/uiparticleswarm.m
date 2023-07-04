function guidata = uiparticleswarm(costFcn,modelspec,init,lb,ub,varargin)
%UIPARTICLESWARM Construct a GUI for running PSO interactively.
%
% -- Usage --
% data = uiparticleswarm(costFcn,modelspec,init,lb,ub)
%
% -- Changelog --
% 2023.07.03 | Created | Wesley Hileman <whileman@uccs.edu>

parser = inputParser;
parser.addRequired('costFcn',@(x)isa(x,'function_handle'));
parser.addRequired('modelspec',@isstruct);
parser.addRequired('init',@isstruct);
parser.addRequired('lb',@isstruct);
parser.addRequired('ub',@isstruct);
parser.addParameter('PopulationSize',500);
parser.addParameter('MaximumIterations',100);
parser.addParameter('PlotInitializeFcn', ...
    @(varargin)[],@(x)isa(x,'function_handle'));
parser.addParameter('PlotUpdateFcn', ...
    @(varargin)[],@(x)isa(x,'function_handle'));
parser.addParameter('WindowName','Particle Swarm Optimization');
parser.parse(costFcn,modelspec,init,lb,ub,varargin{:});
arg = parser.Results;  % structure of validated arguments

costFcn = costFunctionFactory(arg.costFcn,arg.modelspec);
modelspec = arg.modelspec;

guiarg = rmfield(arg,{'costFcn','modelspec','init','lb','ub'});
gui = fastopt.PSOApp(modelspec,arg.init,arg.lb,arg.ub,guiarg);
while true
    % Fetch optimization parameters from user in GUI.
    guidata = gui.run();
    lb = fastopt.pack(guidata.lb,modelspec);
    ub = fastopt.pack(guidata.ub,modelspec);
    init = fastopt.pack(guidata.init,modelspec);
    estimate = init;
    popsize = guidata.PopulationSize;
    maxiter = guidata.MaximumIterations;
    if guidata.stop
        % User canceled the optimization.
        break;
    end
    [plotFcn, stopFcn] = PSOPlotFcnFactory( ...
        arg,guidata);
    gui.lock(stopFcn);  % prevent changes to ui while optimization is running

    % Collect PSO options.
    swarm = zeros(popsize,modelspec.nvars);
    swarm(1,:) = init;
    for k = 1:modelspec.nvars
        swarm(2:end,k) = (ub(k)-lb(k))*rand(popsize-1,1)+lb(k);
    end
    options = optimoptions(@particleswarm,...
        'Display','iter', ...
        'UseParallel',true,...
        'FunctionTolerance',1e-15, ...
        'SwarmSize',popsize,...
        'InitialSwarmMatrix',swarm,...
        'MaxIterations',maxiter, ...
        'MaxStallIterations',200,...
        'FunValCheck','off', ...
        'OutputFcn',plotFcn);  % don't use 'PlotFcn' option, has issues!

    % Run PSO.
    estimate = particleswarm(costFcn,modelspec.nvars,lb,ub,options);

    % Update GUI with new estimate.
    gui.update(fastopt.unpack(estimate,arg.modelspec));
    gui.unlock();  % allow changes to ui now that PSO finished
end

% Collect output data.
guidata = struct;
guidata.values = fastopt.unpack(estimate,modelspec);
guidata.lb = fastopt.unpack(lb,modelspec);
guidata.ub = fastopt.unpack(ub,modelspec);
guidata.arg = arg;
guidata.origin__ = 'fastopt.uiparticleswarm';

end

function [plotFcn, stopFcn] = PSOPlotFcnFactory(arg,guidata)
    plotUpdateIntervalSec = 0.5;
    modelspec = arg.modelspec;
    updatePlotFcn = arg.PlotUpdateFcn;
    updateControlsFcn = guidata.updateControlsFcn;
    updateStatusFcn = guidata.updateStatusFcn;
    updateMark = tic;
    bestJ = Inf;
    stop = false;
    plotFcn = @PSOPlot;
    stopFcn = @PSOStop;

    function halt = PSOPlot(optimValues,state)
        halt = stop;
        if halt
            return;
        end
        if ~strcmp(state,'iter')
            return;
        end
        if toc(updateMark) > plotUpdateIntervalSec
            if optimValues.bestfval < bestJ
                bestJ = optimValues.bestfval;
                bestEstimate = fastopt.unpack(optimValues.bestx,modelspec);
                updatePlotFcn(bestEstimate);
                updateControlsFcn(bestEstimate);
                updateMark = tic;
            end
        end % if
        statusString = sprintf( ...
            "Iteration: %d   BestJ: %.3g   Evals: %d", ...
            optimValues.iteration,optimValues.bestfval,optimValues.funccount);
        updateStatusFcn(statusString);
        drawnow;
    end % PSOPlot()

    function PSOStop()
        % Send stop signal.
        stop = true;
    end
end

function wrapped = costFunctionFactory(costFcn,modelspec)
    wrapped = @costFnWrapper;
    function Jcost = costFnWrapper(vect)
        unpacked = fastopt.unpack(vect,modelspec);
        Jcost = costFcn(unpacked);
    end
end

