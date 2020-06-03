function [segmentInfo, observation, fhmmIntermediate] = segmentModifyZvcMethod_runModel(observation, model, sysparam, fhmmIntermediate)
% use the model to determine the segments of new obs
% [normPos, normVelo] = segmentModiftZvcMethod_preamble(data);
% 
% observation.obsTime = data.time;
% observation.obsData = normPos;
% observation.obsVelo = normVelo;

[segResults, observation, fhmmIntermediate, algSegmentResult] = fgSegment(observation, model, sysparam, fhmmIntermediate);

if ~isempty(segResults)
    
    for i = 1:length(segResults(:, 1))
        ack{i} = '';
    end
    
    segmentInfo.startTimeVal = observation.obsTime(segResults(:, 1));
    segmentInfo.startTimeInd = segResults(:, 1);
    segmentInfo.endTimeVal = observation.obsTime(segResults(:, 1));
    segmentInfo.endTimeInd = segResults(:, 2);
    segmentInfo.segmentName = ack;
    segmentInfo.segmentInclude = ones(size(segResults(:, 1)));
else
    segmentInfo.startTimeVal =[];
    segmentInfo.startTimeInd = [];
    segmentInfo.endTimeVal =[];
    segmentInfo.endTimeInd = [];
    segmentInfo.segmentName =[];
    segmentInfo.segmentInclude = [];
end