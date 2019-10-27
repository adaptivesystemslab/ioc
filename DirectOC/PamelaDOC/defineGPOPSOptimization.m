function setup = defineGPOPSOptimization(config, continuousFunc)

    numDof = length(config.initialState.jointAngles);
    name = sprintf("Robot_%s_DOC", config.name);

    %% Setup for Problem Bounds 
    iphase = 1;
    bounds.phase(iphase).initialtime.lower   =  config.initialState.time;
    bounds.phase(iphase).initialtime.upper   =  config.initialState.time;
    bounds.phase(iphase).finaltime.lower     =  config.finalState.time;
    bounds.phase(iphase).finaltime.upper     =  config.finalState.time;
    bounds.phase(iphase).initialstate.lower  =  [config.initialState.jointAngles', config.initialState.angularVelocities'];
    bounds.phase(iphase).initialstate.upper  =  [config.initialState.jointAngles', config.initialState.angularVelocities'];    
    bounds.phase(iphase).state.lower         =  config.bounds.minState';
    bounds.phase(iphase).state.upper         =  config.bounds.maxState';
    bounds.phase(iphase).finalstate.lower    =  [config.finalState.jointAngles', config.finalState.angularVelocities'];
    bounds.phase(iphase).finalstate.upper    =  [config.finalState.jointAngles', config.finalState.angularVelocities'];
    bounds.phase(iphase).control.lower       =  config.bounds.minControl';
    bounds.phase(iphase).control.upper       =  config.bounds.maxControl';
    bounds.phase(iphase).integral.lower      =  0;
    bounds.phase(iphase).integral.upper      =  1e20;
   
    

    %% Provide Guess of Solution by defining range of values 
%     guess.phase(iphase).state=[config.guess.jointAngles, config.guess.angularVelocities];
%     guess.phase(iphase).control  = config.guess.control;
%     guess.phase(iphase).time     = config.guess.time';
%     guess.phase(iphase).integral = 0;
    guess.phase(iphase).state=3*rand(size([config.guess.jointAngles, config.guess.angularVelocities]));
    guess.phase(iphase).control  = 10*rand(size(config.guess.control));
    guess.phase(iphase).time     = config.guess.time';
    guess.phase(iphase).integral = 0;



    %% Provide Mesh Refinement Method and Initial Mesh 
    mesh.method          = 'hp-LiuRao-Legendre'; %'hp-LiuRao-Legendre';%'hp-PattersonRao';
    mesh.tolerance       = 1e-6;
    mesh.maxiterations   = 10;
    mesh.colpointsmin    = 2;
    mesh.colpointsmax    = 100;
    mesh.phase.colpoints = 100;
    mesh.phase.fraction  = 1;


    %% Assemble Information into Problem Structure 
    setup.mesh                            = mesh;
    setup.name                            = char(name);
    setup.functions.endpoint              = @humanDynEndFunc;
    setup.functions.continuous            = continuousFunc;
    setup.displaylevel                    = 2;
    setup.bounds                          = bounds;
    setup.guess                           = guess;
    setup.nlp.solver                      = 'ipopt';
    setup.nlp.ipoptoptions.linear_solver  = 'ma57';
    setup.nlp.ipoptoptions.tolerance      = 1e-3;
    setup.nlp.ipoptoptions.maxiterations  = 100;
    setup.derivatives.supplier            = 'sparseFD';   % valid options are 'sparseCD', 'sparseFD', 'sparseBD','analytic', and 'adigator'
    setup.derivatives.derivativelevel     = 'first';
    setup.derivatives.dependencies        = 'sparse';
    setup.method                          = 'RPMIntegration';
    setup.scales.method                   = 'automatic-bounds';
end
