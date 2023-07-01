function [model, trajectory] = particleswarm(costFn,modelspec,lb,ub,varargin)
    %PARTICLESWARM

    parser = inputParser;
    parser.addParameter('display','iter');
    parser.addParameter('particleCount',1000);
    parser.addParameter('swarmIterations',200);
    parser.addParameter('fminconIterations',1000);
    parser.addParameter('fminconFunEvals',1e9);
    parser.addParameter('hybrid','fmincon');
    parser.addParameter('trackTrajectory',false);
    parser.addParameter('initial',[]);
    parser.parse(varargin{:});

    display = parser.Results.display;
    pcntSwarm = parser.Results.particleCount;
    iterSwarm = parser.Results.swarmIterations;
    iterFmincon = parser.Results.fminconIterations;
    funEvalsFmincon = parser.Results.fminconFunEvals;
    hybrid = parser.Results.hybrid;
    trackTrajectory = parser.Results.trackTrajectory;
    initial = parser.Results.initial;

    % Prepare particle swarm matrix. Each row specifies the 
    % parameters associated with a particle in the swarm.
    swarm = zeros(pcntSwarm,modelspec.nvars);
    if ~isempty(initial)
        swarm(1,:) = initial;
        for k = 1:modelspec.nvars
            swarm(2:end,k) = (ub(k)-lb(k))*rand(pcntSwarm-1,1)+lb(k);
        end
    else
        for k = 1:modelspec.nvars
            swarm(:,k) = (ub(k)-lb(k))*rand(pcntSwarm,1)+lb(k);
        end
    end
    
    if strcmpi(hybrid,'fmincon')
        optionsHybrid = optimoptions(@fmincon,'Display',display,...
            'MaxFunEvals',funEvalsFmincon,'MaxIter',iterFmincon,...
            'TolFun',1e-15,'TolX',1e-15,'TolCon',1e-15, ...
            'OutputFcn',@iterfcnFmincon);
        options = optimoptions(@particleswarm,...
            'Display',display,'UseParallel',true,...
            'FunctionTolerance',1e-15,'SwarmSize',pcntSwarm,...
            'InitialSwarmMatrix',swarm,...
            'MaxIterations',iterSwarm,'MaxStallIterations',200,...
            'FunValCheck','off','OutputFcn',@iterfcnParticleSwarm,...
            'HybridFcn',{@fmincon,optionsHybrid});
        maxIterations = iterSwarm+iterFmincon;
    else
        options = optimoptions(@particleswarm,...
            'Display',display,'UseParallel',true,...
            'FunctionTolerance',1e-15,'SwarmSize',pcntSwarm,...
            'InitialSwarmMatrix',swarm,...
            'MaxIterations',iterSwarm,'MaxStallIterations',200,...
            'FunValCheck','off','OutputFcn',@iterfcnParticleSwarm);
        maxIterations = iterSwarm;
    end
    
    % Set up storage for tracking parameter trajectory during the
    % optimization.
    if trackTrajectory
        paramvalues = zeros(modelspec.nvars,maxIterations);
        step = zeros(1,maxIterations);
    else
        paramvalues = [];
        step = [];
    end
    cursor = 1;

    % Run optimization.
    vect = particleswarm(@costFnWrapper,modelspec.nvars,lb,ub,options);
    model = fastopt.unpack(vect,modelspec);

    if trackTrajectory
        % Convert matrix of column vectors to structure array of models.
        paramvalues = cellfun( ...
            @(x)fastopt.unpack(x,modelspec,'sparse',true,'flat',true), ...
            num2cell(paramvalues(:,1:cursor-1),1));
        trajectory.paramvalues = paramvalues;
        trajectory.step = step(1:cursor-1);
        trajectory.iteration = 1:cursor-1;
    end

    function Jcost = costFnWrapper(vect)
        unpacked = fastopt.unpack(vect,modelspec);
        Jcost = costFn(unpacked);
    end

    function stop = iterfcnParticleSwarm(optimOptions,state)
        stop = false;
        if strcmpi(state,'iter')
            if trackTrajectory
                paramvalues(:,cursor) = optimOptions.bestx;
                step(cursor) = 1;
                cursor = cursor + 1;
            end
        end
    end

    function stop = iterfcnFmincon(x,~,state)
        stop = false;
        if strcmpi(state,'iter')
            if trackTrajectory
                paramvalues(:,cursor) = x;
                step(cursor) = 2;
                cursor = cursor + 1;
            end
        end
    end

end