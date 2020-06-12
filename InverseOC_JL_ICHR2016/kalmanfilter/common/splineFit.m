classdef splineFit < handle
    properties
        splineViaPoints = [];
        
        splineDof = [];
        splineTimeLength = [];
        
        splineType = [];
        splineParam = []; % 
    end
    
    methods
        function obj = splineFit(splineType)
            obj.splineType = splineType;
        end
        
        function obj = setViaPointsStruct(obj, splineViaPoints)
            % set the struct version
            obj.splineViaPoints = splineViaPoints;

            % proprogate the chage to the vector version
            obj.splineParam = [];
        end
        
        function obj = updateViaPoints(obj, q_mat)
            obj.splineViaPoints.q = q_mat{1};
            obj.splineViaPoints.dq = q_mat{2};
            obj.splineViaPoints.ddq = q_mat{3};
        end
        
        function q_mat = downloadViaPoints(obj)
            q_mat{1} = obj.splineViaPoints.q;
            q_mat{2} = obj.splineViaPoints.dq;
            q_mat{3} = obj.splineViaPoints.ddq;
        end
        
%         function obj = setViaPointsCell(splineViaPointsCell)
%             % set the struct version
%             obj.splineViaPoints = splineViaPoints;
%             
%             % proprogate the chage to the vector version
%             obj.splineParam = [];
%         end

        function newSplineViaPoints = splitSplineViaPoints(obj, ind1, ind2)
            newSplineViaPoints.inds = obj.splineViaPoints.inds(ind1:ind2);
            newSplineViaPoints.time = obj.splineViaPoints.time(ind1:ind2);
            newSplineViaPoints.q = obj.splineViaPoints.q(ind1:ind2, :);
            newSplineViaPoints.dq = obj.splineViaPoints.dq(ind1:ind2, :);
            newSplineViaPoints.ddq = obj.splineViaPoints.ddq(ind1:ind2, :);
        end

        function obj = calcSplineParam(obj)
            % pull out the spline points:
            %  t: time value corresponding to the key points
            %  y: magnitude value corresponding to the key points
            switch obj.splineType
                case enumSplineTypes.piecewiseQuintic
                    for i = 1:length(obj.splineViaPoints.inds)-1
                        newSpline = splineFit(enumSplineTypes.quinticSpline);
                        newSpline.setViaPointsStruct(obj.splitSplineViaPoints(i, i+1));
                        sp_qd0{i} = newSpline.calcSplineParam();
                        sp_qd1{i} = [];
                        sp_qd2{i} = [];
                        sp_qd3{i} = [];
                    end
                    
                    time = [];
                    inds = [];
                    
                case enumSplineTypes.quinticSpline
                    for i = 1:size(obj.splineViaPoints(1).q, 2)
                        t0 = obj.splineViaPoints.time(1);
                        q0 = obj.splineViaPoints.q(1, i);
                        dq0 = obj.splineViaPoints.dq(1, i);
                        ddq0 = obj.splineViaPoints.ddq(1, i); 
                        
                        tf = obj.splineViaPoints.time(end); 
                        qf = obj.splineViaPoints.q(end, i); 
                        dqf = obj.splineViaPoints.dq(end, i); 
                        ddqf = obj.splineViaPoints.ddq(end, i);
                        
                        sp_qd0{i} = quinticPoly(q0, dq0, ddq0, qf, dqf, ddqf, t0, tf);
                        sp_qd1{i} = polyder_num(sp_qd0{i});
                        sp_qd2{i} = polyder_num(sp_qd1{i});
                        sp_qd3{i} = polyder_num(sp_qd2{i});
                    end
                    
                    time = obj.splineViaPoints.time;
                    inds = obj.splineViaPoints.inds;
            end
            
            obj.splineParam.t = time;
            obj.splineParam.inds = inds;
            obj.splineParam.sp_qd0 = sp_qd0;
            obj.splineParam.sp_qd1 = sp_qd1;
            obj.splineParam.sp_qd2 = sp_qd2;
            obj.splineParam.sp_qd3 = sp_qd3;
        end
        
        function [q, dq, ddq] = calcQ(obj, time)
            switch obj.splineType
                case enumSplineTypes.piecewiseQuintic
                    dofCount = length(obj.splineParam.sp_qd0{1}.splineParam.sp_qd0);
                    tCount = length(time);
                    
                    q = zeros(tCount, dofCount);
                    dq = zeros(tCount, dofCount);
                    ddq = zeros(tCount, dofCount);
                    
                    for ind_n = 1:length(obj.splineParam.sp_qd0) 
                        [q_temp, dq_temp, ddq_temp] = obj.splineParam.sp_qd0{ind_n}.calcQ(time);
                        
                        ind_start = obj.splineParam.sp_qd0{ind_n}.splineParam.inds(1);
                        ind_end = obj.splineParam.sp_qd0{ind_n}.splineParam.inds(end);
                        
                        switch ind_n
                            case 1
                                
                            otherwise
                                ind_start = ind_start + 1; % this keeps it causal
                        end
                        
                        % the quintic is continuous, so need to pull out the
                        % relavent section of the spline
                        q(ind_start:ind_end, :) = q_temp(ind_start:ind_end, :);
                        dq(ind_start:ind_end, :) = dq_temp(ind_start:ind_end, :);
                        ddq(ind_start:ind_end, :) = ddq_temp(ind_start:ind_end, :);
                    end

                case enumSplineTypes.quinticSpline
                    dofCount = length(obj.splineParam.sp_qd0);
                    tCount = length(time);
                    
                    q = zeros(tCount, dofCount);
                    dq = zeros(tCount, dofCount);
                    ddq = zeros(tCount, dofCount);
                    
                    for ind_n = 1:length(obj.splineParam.sp_qd0)
                        q(:, ind_n) =     polyval(obj.splineParam.sp_qd0{ind_n}, time)';
                        dq(:, ind_n) =    polyval(obj.splineParam.sp_qd1{ind_n}, time)';
                        ddq(:, ind_n) =   polyval(obj.splineParam.sp_qd2{ind_n}, time)';
                    end
            end
        end
    end
end



