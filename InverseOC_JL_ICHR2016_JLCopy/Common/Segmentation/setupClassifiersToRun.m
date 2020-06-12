function classifierSetToRun = setupClassifiersToRun(classifierPackage)
% set up the different classifier combinations that we are testing for

counter = 0;
classiferIterationCount = length(classifierPackage.aggregatorList)*length(classifierPackage.dimReductionList)*length(classifierPackage.classifierList);
classifierSetToRun = cell(classiferIterationCount, 3);
% assemble the classifier combination
for ind_dimReduct = 1:length(classifierPackage.dimReductionList)
    for ind_aggregator = 1:length(classifierPackage.aggregatorList)
        for ind_classifier = 1:length(classifierPackage.classifierList)
            counter = counter + 1;
            classifierSetToRun{counter, 1} = classifierPackage.dimReductionList{ind_dimReduct};
            classifierSetToRun{counter, 2} = classifierPackage.classifierList{ind_classifier};
            classifierSetToRun{counter, 3} = classifierPackage.aggregatorList{ind_aggregator};
        end
    end
end