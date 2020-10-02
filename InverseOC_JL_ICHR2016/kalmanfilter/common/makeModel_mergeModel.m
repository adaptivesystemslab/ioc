function [modelStructTrc, modelStructImu] = makeModel_mergeModel(overallModelStruct)
    % copy out the link offsets for kinematic transform 
    kinT = {};
    for i = 1:length(overallModelStruct)
        for j = 1:length(overallModelStruct(1).kinematicTransform)
            kinT{j}{i} = overallModelStruct(i).kinematicTransform(j).t;
        end
    end
    for j = 1:length(overallModelStruct(1).kinematicTransform)
        kinTransOut(j).frameName = overallModelStruct(1).kinematicTransform(j).frameName;
        kinTransOut(j).t = meanSVD_T(kinT{j});
    end
    
    dynM = {};
    dynCom = {};
    dynI = {};
    for i = 1:length(overallModelStruct)
        for j = 1:length(overallModelStruct(1).dynamicTransform)
            dynM{j}{i} = overallModelStruct(i).dynamicTransform(j).m;
            dynCom{j}{i} = overallModelStruct(i).dynamicTransform(j).com;
            dynI{j}{i} = overallModelStruct(i).dynamicTransform(j).I;
        end
    end
    for j = 1:length(overallModelStruct(1).dynamicTransform)
        dynTransOut(j).frameName = overallModelStruct(1).dynamicTransform(j).frameName;
        dynTransOut(j).m = meanR(dynM{j});
        dynTransOut(j).com = meanR(dynCom{j});
        dynTransOut(j).I = meanR(dynI{j});
    end
    
    senT = {};
    for i = 1:length(overallModelStruct)
        for j = 1:length(overallModelStruct(1).sensorTransform)
            senT{j}{i} = overallModelStruct(i).sensorTransform(j).t;
        end
    end
    for j = 1:length(overallModelStruct(1).sensorTransform)
        sensTransOut(j).frameName = overallModelStruct(1).sensorTransform(j).frameName;
        sensTransOut(j).sensorName = overallModelStruct(1).sensorTransform(j).sensorName;
        sensTransOut(j).decorators = overallModelStruct(1).sensorTransform(j).decorators;
        
        if isempty(sensTransOut(j).frameName)
            sensTransOut(j).frameName = overallModelStruct(2).sensorTransform(j).frameName;
            sensTransOut(j).sensorName = overallModelStruct(2).sensorTransform(j).sensorName;
            sensTransOut(j).decorators = overallModelStruct(2).sensorTransform(j).decorators;
        end
        
        sensTransOut(j).t = meanSVD_T(senT{j});
    end

    senSecT = {};
    for i = 1:length(overallModelStruct)
        for j = 1:length(overallModelStruct(1).sensorSecondaryTransformTrc)
            senSecT{j}{i} = overallModelStruct(i).sensorSecondaryTransformTrc(j).t;
        end
    end
    for j = 1:length(overallModelStruct(1).sensorSecondaryTransformTrc)
        sensSecTransOutTrc(j).frameName = overallModelStruct(1).sensorSecondaryTransformTrc(j).frameName;
        sensSecTransOutTrc(j).sensorName = overallModelStruct(1).sensorSecondaryTransformTrc(j).sensorName;
        sensSecTransOutTrc(j).decorators = overallModelStruct(1).sensorSecondaryTransformTrc(j).decorators;
        
        if isempty(sensSecTransOutTrc(j).frameName)
            sensSecTransOutTrc(j).frameName = overallModelStruct(2).sensorSecondaryTransformTrc(j).frameName;
            sensSecTransOutTrc(j).sensorName = overallModelStruct(2).sensorSecondaryTransformTrc(j).sensorName;
            sensSecTransOutTrc(j).decorators = overallModelStruct(2).sensorSecondaryTransformTrc(j).decorators;
        end
        
        sensSecTransOutTrc(j).t = meanSVD_T(senSecT{j});
    end
    
    modelStructTrc.kinematicTransform = kinTransOut;
    modelStructTrc.dynamicTransform = dynTransOut;
    modelStructTrc.sensorTransform = sensTransOut;
    modelStructTrc.sensorSecondaryTransform = sensSecTransOutTrc;
    
    if isfield(overallModelStruct(1), 'sensorSecondaryTransformImu') && ~isempty(overallModelStruct(1).sensorSecondaryTransformImu)
        senSecT = {};
        for i = 1:length(overallModelStruct)
            for j = 1:length(overallModelStruct(1).sensorSecondaryTransformImu)
                senSecT{j}{i} = overallModelStruct(i).sensorSecondaryTransformImu(j).t;
            end
        end
        for j = 1:length(overallModelStruct(1).sensorSecondaryTransformImu)
            sensSecTransOutImu(j).frameName = overallModelStruct(1).sensorSecondaryTransformImu(j).frameName;
            sensSecTransOutImu(j).sensorName = overallModelStruct(1).sensorSecondaryTransformImu(j).sensorName;
            sensSecTransOutImu(j).decorators = overallModelStruct(1).sensorSecondaryTransformImu(j).decorators;
            
            if isempty(sensSecTransOutImu(j).frameName)
                sensSecTransOutImu(j).frameName = overallModelStruct(2).sensorSecondaryTransformImu(j).frameName;
                sensSecTransOutImu(j).sensorName = overallModelStruct(2).sensorSecondaryTransformImu(j).sensorName;
                sensSecTransOutImu(j).decorators = overallModelStruct(2).sensorSecondaryTransformImu(j).decorators;
            end
            
            [sensSecTransOutImu(j).t, sensSecTransOutImu(j).dist] = meanSVD_T(senSecT{j});
        end
        
        modelStructImu.kinematicTransform = kinTransOut;
        modelStructImu.dynamicTransform = dynTransOut;
        modelStructImu.sensorTransform = sensTransOut;
        modelStructImu.sensorSecondaryTransform = sensSecTransOutImu;
    else
        modelStructImu = [];
    end
end