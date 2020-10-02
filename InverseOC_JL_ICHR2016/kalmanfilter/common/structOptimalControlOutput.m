classdef structOptimalControlOutput < handle
    properties
        time;
        inds;
        q;
        c;
        ocStruct;
        optStruct;
        
        optPivot;
        
        splinePoints;
        constPoints;
        resnorm;
    end
    
    methods
         function obj = structOptimalControlOutput()
             
         end
        
        function obj = loadData(varargin)
            obj.time = varargin{1};
            obj.q = varargin{2};
            obj.c = varargin{3};
            obj.optStruct = varargin{4};
            obj.resnorm = varargin{5};
        end
        
        function obj = saveFeatures(obj, featureWin, currWin)
            obj.time = featureWin.time;
            obj.q = featureWin.q;
            obj.inds = currWin;
        end
        
        function obj = saveOCStruct(obj, ocStruct)
            obj.ocStruct = ocStruct;
            obj.optStruct = ocStruct.getIocOptStruct();
            
            obj.c = ocStruct.getIocOptStructBestPivot().c;
            obj.resnorm = ocStruct.getIocOptStructBestPivot().resnorm;
            
            obj.optPivot = obj.ocStruct.costFunctionStruct(obj.ocStruct.optPivotInd).name;
            
            obj.splinePoints = obj.inds(ocStruct.featureSet.splineModel.splineViaPoints.inds);
            obj.constPoints = obj.inds(ocStruct.constPoints.inds);
        end
        
        function obj = writeToFile(obj, savePath)
            writeStruct.runTime = datestr(now);
            writeStruct.t_start = obj.inds(1);
            writeStruct.t_end = obj.inds(end);
            writeStruct.qCk_Rmse = 0;
            writeStruct.qDoc_Rmse = 0;
            
            writeStruct.optpivot = obj.optPivot;
            
            for i = 1:length(obj.ocStruct.costFunctionStruct)
                writeStruct.(['C_gen_cost_' obj.ocStruct.costFunctionStruct(i).name]) = 0;
            end
            
            for i = 1:length(obj.ocStruct.costFunctionStruct)
                writeStruct.(['C_rec_cost_' obj.ocStruct.costFunctionStruct(i).name]) = obj.c(i);
            end
            
            writeStruct.rmse = 0;
            writeStruct.resnorm	= obj.resnorm;

            newTable = struct2table(writeStruct);
            
            % if file exist already, load it 
            if exist(savePath, 'file')
                oldTable = readtable(savePath);
                newTable = [oldTable; newTable];
            end
            
            writetable(newTable, savePath);
        end
    end
end