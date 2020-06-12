function [dimReduct] = dimReductProcessing(dimReductSelect, data)

    switch dimReductSelect.name
        case 'None'
            dimReduct = NoTransform;

            %         case 'ksPCA'
            %             dimReduct = KSPCATransform;
            %
            %             if isfield(dimReductSelect, 'degree')
            %                 dimReduct.init(dimReductSelect.y_kernel, dimReductSelect.x_kernel, ...
            %                     dimReductSelect.sigma, dimReductSelect.degree);
            %             end

        case 'ksPCA'
            dimReduct = KSPCATransform2;

            if isfield(dimReductSelect, 'degree')
                dimReduct.init(dimReductSelect.y_kernel, dimReductSelect.x_kernel, ...
                    dimReductSelect.sigma, dimReductSelect.degree, dimReductSelect.LMatrixMultiplier);
            end

        case 'lasso'
            dimReduct = LassoTransform(dimReductSelect.degree);

        case 'PCA'
            if ~isfield(dimReductSelect, 'degree')
                % generate default settings
                dimReductSelect.degree = 7;
            end

            dimReduct = PCATransform(dimReductSelect.degree);

        case 'fPCA'
            dimReduct = fPCATransform(dimReductSelect.degree);

        case 'FDA'
            if ~isfield(dimReductSelect, 'degree')
                % generate default settings
                dimReductSelect.degree = 2;
            end

            dimReduct = FDATransform(dimReductSelect.degree);

        case 'LMNN'
            dimReduct = LMNNTransform(dimReductSelect.degree);

        case 'LMNN2'
            dimReduct = LMNNTransform2(dimReductSelect.degree);
    end
    
    if isa(dimReduct, 'fPCATransform')
        dimReduct.train(data);
    else
        dimReduct.train(data.trainingData, data.trainingLabel);
    end
    
%     figure
%     hold on
%     plot(testingData(1, :), testingData(2, :), '*k')
%     find1 = find(testingLabel == 1);
%     plot(testingData(1, find1), testingData(2, find1), '*r')
%     hold on
%     find0 = find(testingLabel == 0);
%     plot(testingData(1, find0), testingData(2, find0), '*b')
%     grid on
end