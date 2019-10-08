function setup = defineGPOPSOptimization(config, continuousFunction)

    numDof = length(config.initialState.jointAngles);
    name = sprintf("Robot_%s_DOC", config.name);

    %% Setup for Problem Bounds 
    iphase = 1;
    bounds.phase(iphase).initialtime.lower   =  config.initialState.time;
    bounds.phase(iphase).initialtime.upper   =  config.initialState.time;
    bounds.phase(iphase).finaltime.lower     =  config.finalState.time;
    bounds.phase(iphase).finaltime.upper     =  config.finalState.time;
    bounds.phase(iphase).initialstate.lower  =  [arrayfun(@(x) deg2rad(x), config.initialState.jointAngles)', config.initialState.angularVelocities'];
    bounds.phase(iphase).initialstate.upper  =  [arrayfun(@(x) deg2rad(x), config.initialState.jointAngles)', config.initialState.angularVelocities'];    
    bounds.phase(iphase).state.lower         =  [arrayfun(@(x) deg2rad(x), config.bounds.minState(1:numDof))', config.bounds.minState(numDof+1:end)'];
    bounds.phase(iphase).state.upper         =  [arrayfun(@(x) deg2rad(x), config.bounds.maxState(1:numDof))', config.bounds.maxState(numDof+1:end)'];
    bounds.phase(iphase).finalstate.lower    =  [arrayfun(@(x) deg2rad(x), config.finalState.jointAngles)', config.finalState.angularVelocities'];
    bounds.phase(iphase).finalstate.upper    =  [arrayfun(@(x) deg2rad(x), config.finalState.jointAngles)', config.finalState.angularVelocities'];
    bounds.phase(iphase).control.lower       =  config.bounds.minControl';
    bounds.phase(iphase).control.upper       =  config.bounds.maxControl';
    bounds.phase(iphase).integral.lower      =  0;
    bounds.phase(iphase).integral.upper      =  1;
    
    auxdata.state.lower = bounds.phase(iphase).state.lower;
    auxdata.state.upper = bounds.phase(iphase).state.upper;
    auxdata.control.lower = bounds.phase(iphase).control.lower;
    auxdata.control.upper = bounds.phase(iphase).control.upper;
    auxdata.integral.lower = bounds.phase(iphase).integral.lower;
    auxdata.integral.upper = bounds.phase(iphase).integral.upper;
    

    %% Provide Guess of Solution by defining range of values
    anglesGuess = []; velocityGuess = []; controlGuess = [];
    for i=1:length(config.guess)
        temp = config.guess(i);
        anglesGuess = [anglesGuess, arrayfun(@(x) deg2rad(x), temp.jointAngles)];
        velocityGuess = [velocityGuess, temp.angularVelocities];
        controlGuess = [controlGuess, temp.control];    
    end
    
    if isfield(config.guess(1), 'time')
        timeGuess = config.guess(1).time';
    else
        timeGuess = [config.initialState.time; config.finalState.time];
    end
    
    guess.phase(iphase).state    = [anglesGuess, velocityGuess];
    guess.phase(iphase).control  = controlGuess;
    guess.phase(iphase).time     = timeGuess;
    guess.phase(iphase).integral = 0;


    %% Provide Mesh Refinement Method and Initial Mesh 
    mesh.method          = 'hp-LiuRao-Legendre';%'hp-PattersonRao';
    mesh.tolerance       = 1e-3;
    mesh.maxiterations   = 5;
    mesh.colpointsmin    = 2;
    mesh.colpointsmax    = 10;
    mesh.phase.colpoints = 2*ones(1,5);
    mesh.phase.fraction  = ones(1,5)/5;

    %% Assemble Information into Problem Structure 
    setup.mesh                            = mesh;
    setup.name                            = char(name);
    setup.functions.endpoint              = @robotArmEndpoint;
    setup.functions.continuous            = continuousFunction;
    setup.displaylevel                    = 2;
    setup.bounds                          = bounds;
    setup.guess                           = guess;
    setup.auxdata                         = auxdata;
    setup.nlp.solver                      = 'ipopt';
    setup.nlp.ipoptoptions.linear_solver  = 'ma57';
    setup.nlp.ipoptoptions.tolerance      = 1e-3;
    setup.nlp.ipoptoptions.maxiterations  = 25;
    setup.derivatives.supplier            = 'sparseCD';   % valid options are 'sparseCD', 'sparseFD', 'sparseBD','analytic', and 'adigator'
    setup.derivatives.derivativelevel     = 'second';
    setup.derivatives.dependencies        = 'sparse';
    setup.method                          = 'RPM-Differentiation';
    setup.scales.method                   = 'automatic-bounds';
end

