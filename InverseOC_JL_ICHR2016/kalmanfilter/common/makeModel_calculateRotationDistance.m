function [modelStructTrc, modelStructImu] = makeModel_calculateRotationDistance(overallModelStruct)
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
        
        sensSecTransOutImu(j).t = meanT(senSecT{j});
    end
    
    modelStructTrc.kinematicTransform = kinTransOut;
    modelStructTrc.dynamicTransform = dynTransOut;
    modelStructTrc.sensorTransform = sensTransOut;
    modelStructTrc.sensorSecondaryTransform = sensSecTransOutTrc;
    
    if length(overallModelStruct(1).sensorSecondaryTransformImu) > 0
        modelStructImu.kinematicTransform = kinTransOut;
        modelStructImu.dynamicTransform = dynTransOut;
        modelStructImu.sensorTransform = sensTransOut;
        modelStructImu.sensorSecondaryTransform = sensSecTransOutImu;
    else
        modelStructImu = [];
    end
end