classdef ModelJoint < handle
    %ARTICULATEDJOINT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Kinematic parameters
        twist = 0; % alpha
        offset = 0; % d
        length = 0; % a/r
        angleOffset = 0; % theta
        endEffector = 0;
        
        % State variables
        q = 0; dq = 0; ddq = 0;
        
        % Dynamic parameters
        mass = 0;
        com = [0 0 0] % center of mass
        inertiaTensor = zeros(3,3);
    end
    
    methods
                
        function obj = ModelJoint(params)
            obj.twist = params.a;
            obj.offset = params.d;
            obj.length = params.r;
            obj.angleOffset = params.theta;
            obj.endEffector = params.endEffector;
                            
            % If a physical link is associated to this joint
            if params.virtual == 0
                obj.mass = params.m;            
                obj.com = params.com;
                obj.inertiaTensor = computeInertiaTensor(obj.mass,...
                    obj.com, params.tensor);                
            end
        end       
        
        function updateState(obj, newQ, newDQ)
            obj.q = newQ;
            obj.dq = newDQ;
        end
        
    end
end

