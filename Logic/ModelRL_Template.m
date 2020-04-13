classdef ModelRL_Template < ModelRL
    properties
    end
    
    methods
        function loadModel(obj, xmlPath, trialInfo)
            % TO IMPLEMENT
            % This function implements the obj.model object. To do this,
            % an XML-based transformation file is to be constructed and
            % loaded via:
            %     model = rlCModel(rlModelPath);
            % If required, the link length and inertial parameters should
            % be populated here as well
            
            % EXAMPLE CODE (or consult ModelRL_Jumping_KW and other
            % ModelRLs for examples) to populate a single joint of a full
            % body model
            
            % load the XML file to create the RL model. 
            model = rlCModel(xmlPath); % load the XML file
            
            % define some parameters to apply to the RL model. For this
            % example, we will apply the transform to the upper right leg
            currTransformName = 'length_rhip_rknee'; % a transformation tag tying two frames together. there is a corresponding entry in the .xml
            currBodyName = 'body_rhip_rknee'; % a body tag that can hold mass/inertial param definitions. there is a corresponding entry in the .xml
            dumasFrameStr = 'thigh'; % the Dumas name for this frame
            linkLength = 0.50; 
            gender = 'M';
            weight = 60;
            
            % full array of names, used to locate the proper frame to apply
            % these modifications
            allTransformNames = {model.transforms.name};
            allBodyNames = {model.bodies.name};
            
            % now apply the transform to the link length
            trNum = find(ismember(allTransformNames,currTransformName)==1); % find the index of the transform we want to modify
            T = model.transforms(trNum).t;
            FT = model.transforms(trNum).frame_in.t;
            T(1:3,4) = FT(1:3,1:3)'*[0,0,linkLength]'; % apply the transform
            model.transforms(trNum).t = T; % save the updated transform
            model.forwardPosition; % proprogate the change down the kinematic chain
            
            rotMtxDumas = rotx(pi/2); % a rotation matrix to get the Dumas frame into our local mocap frame
            mass = lookupTableDumas('mass', dumasFrameStr, gender, [], [])*double(weight); % use Dumas2007, an anthropometric paper, to get mass/com/inertia
            comScale = lookupTableDumas('com', dumasFrameStr, gender, norm(linkLength), []);
            inertialScale = lookupTableDumas('inertial', dumasFrameStr, gender, norm(linkLength), mass);
            com = rotMtxDumas*comScale;
            inertial = rotMtxDumas*inertialScale;
            bodyNum = find(ismember(allBodyNames, currBodyName) == 1);
            model.bodies(bodyNum).m = mass;
            model.bodies(bodyNum).com = com;
            model.bodies(bodyNum).I = inertial;
            
            obj.model = model;
            obj.model.forwardPosition();
        end
        
        function traj = loadData(obj, trialInfo)
            % TO IMPLEMENT
            % this function loads trajectory data to be used to calculate
            % the IOC. The format expected is as follows:
            
            % traj.q (n x d) - joint angle data, in the order that RL
            %   expects it (ie {model.joints.name})
            % traj.dq (n x d) - joint velocity data. can diff(q) to get it 
            % traj.tau (n x d) - joint torque data. can obtain via RL if
            %   have q, dq, and ddq
            % traj.trajT (n x d) - time trajectory data
            % traj.trajX (n x d) - state trajectory data: q and dq
            % traj.trajU (n x d) - control trajectory data: tau
            % traj.frameInds (1 x n) - frames to run for IOC
            
            % EXAMPLE CODE
            dt = 0.01;
            t = 0:dt:8*pi;
            qIndiv = sin(t);
            q = zeros(length(qIndiv), length(model.joints));
            q(:, 7) = qIndiv;
            
            % don't have dq or tau. filtData will generate them accordingly
            traj = obj.filtData(t, dt, q, trialInfo);
        end
        
        function featureSpecialize(obj)
            % if there's any customized specialization that is required,
            % the individual modelRLs should override this
        end
    end
end