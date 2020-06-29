function dataUse = clusterApplyDataFormat(obj, input, output, trainingInputDR, sigdofcellarray)

if ~exist('output', 'var')
    output = [];
end

if ~exist('trainingInputDR', 'var') 
    trainingInputDR = [];
end

if ~exist('sigdofcellarray', 'var')
    sigdofcellarray = [];
end

% each inputuse cell is a 'plate' of data that we will use
[inputUse, currIndUse, dof] = sortDataPlates(obj, input, sigdofcellarray);
        
for ind_plate = 1:length(inputUse)
    % remove the mean from all dim
    if ~isempty(trainingInputDR)
        % application phase
        currPlateMean = trainingInputDR{ind_plate}.mean;
    else
        % training phase
        currPlateMean = removeMeanFct(obj, inputUse{ind_plate});
    end
    currPlate = inputUse{ind_plate} - repmat(currPlateMean, size(inputUse{ind_plate}, 1), 1);
    
    % now for variance
    if ~isempty(trainingInputDR)
        % application phase
        currPlateNorm = trainingInputDR{ind_plate}.normMultiplier;
    else
        % training phase
        currPlateNorm = removeVarianceFct(obj, inputUse{ind_plate});
    end
    currPlate = currPlate ./ repmat(currPlateNorm, size(inputUse{ind_plate}, 1), 1);
    
    switch obj.dataFormatfitting
        case 'AI'
            % leave as is
            inputPlate{ind_plate} = currPlate;
            
        case 'CR'
            % least squares
%             currPlate = currPlate / currPlateNorm;
            inputPlate{ind_plate} = mapToCircle(currPlate, output, 0); % xp = Ppca' * (x-mu) * Pc;
            
%             dataToFitToCircle = [currPlate; currPlate(:, 1) -currPlate(:, 2)];
%             h_la = figure; scatter(dataToFitToCircle(:, 1), dataToFitToCircle(:, 2)); hold on
%             lala = fit_ellipse( dataToFitToCircle(:, 1),dataToFitToCircle(:, 2), h_la );
            
        case 'FI'
            
    end
    
    switch obj.dataFormatQuad
        case 1
            % we only want the data to occupy one quad
            inputPlate{ind_plate} = abs(inputPlate{ind_plate});
         
            
        case 4
            % leave
    end
    
    % if the DOF is 2 or less, we will use the joint angle as the
    % threshold
%     if size(inputPlate{ind_plate}, 2) <= 2

    switch obj.plateSelection
        case 'thresangle'
            [polarTheta, polarRho] = cart2pol(inputPlate{ind_plate}(:, 1), inputPlate{ind_plate}(:, 2));
            dataUse{ind_plate}.data = polarTheta;
            dataUse{ind_plate}.thres = polarRho;
            
        case 'thresmag'
            mags = normVector(inputPlate{ind_plate});
            dataUse{ind_plate}.data = inputPlate{ind_plate};
            dataUse{ind_plate}.thres = mags;
            
        otherwise
            dataUse{ind_plate}.data = inputPlate{ind_plate};
            dataUse{ind_plate}.thres = ones(size(inputPlate{ind_plate}, 1), 1);
    end
    
    dataUse{ind_plate}.dof = dof{ind_plate};
    dataUse{ind_plate}.dofInd = currIndUse{ind_plate};
    dataUse{ind_plate}.mean = currPlateMean;
    dataUse{ind_plate}.normMultiplier = currPlateNorm;
end

if 0
    % circle plots
    h1 = figure('position', [       1384         147        1240         913]);
    
    upperBound = length(inputPlate);
    if upperBound > 3
        upperBound = 3;
    end
    
    for ind = 1:upperBound
        h1(ind) = subplot(1, upperBound, ind);
        
        if ~isempty(dataUse{ind}.dofInd)
            currOutput = output(dataUse{ind}.dofInd);
        else
            currOutput = output;
        end
        
        plot(inputPlate{ind}(currOutput == 1, 1), inputPlate{ind}(currOutput == 1, 2), 'r.'); hold on
        plot(inputPlate{ind}(currOutput == 0, 1), inputPlate{ind}(currOutput == 0, 2), 'b.');
        xlim([-3 3]); ylim([-3 3]);
        linkaxes(h1, 'xy');
    end
end

if 0
    % time series plots
    h1 = figure('position', [       1384         147        1240         913]);
    
    for ind = 1:3
        h1(ind) = subplot(1, upperBound, ind);
        
        if ~isempty(dataUse{ind}.dofInd)
            currOutput = output(dataUse{ind}.dofInd);
        else
            currOutput = output;
        end
        
        t = 1:size(inputPlate{ind}, 1);
        
        plot(t(currOutput == 1), polarTheta(currOutput == 1), 'r.'); hold on
        plot(t(currOutput == 0), polarTheta(currOutput == 0), 'b.');
        %         xlim([-1.5 1.5]); ylim([-1.5 1.5]);
        %         linkaxes(h1, 'xy');
    end
end

end

function currPlateMean = removeMeanFct(obj, currInputUse)
    switch obj.removeMean
        case 'none'
            currPlateMean = zeros(1, size(currInputUse, 2));

        case 'mean'
            % offset is calculated from taking the mean
            currPlateMean = mean(currInputUse);

            switch obj.dataFormatName
                case 'PC'

                case 'PP'
                    switch obj.dataFormatDofStr
                        case {'all', 'alldof', 'sigdof'}
                            % these are expected to be circles, and thus mean
                            % on both dofs is suitable

                        case {'normalldof', 'normsigdof'}
                            % these signals have been normalized, and thus form
                            % half circles
                            currPlateMean(2) = 0;
                    end
            end

        case 'circle'
            % offset is calculated from taking the circle fit
            switch obj.dataFormatName
                case 'PC'
                    dataToFitToCircle = currInputUse;
                    
                case 'PP'
                    switch obj.dataFormatDofStr
                        case {'all', 'alldof', 'sigdof'}
                            dataToFitToCircle = currInputUse;

                        case {'normalldof', 'normsigdof'}
                            % offset is calculated from taking the circle fit
                            dataToUse = currInputUse;
                            dataToFitToCircle = [dataToUse; dataToUse(:, 1) -dataToUse(:, 2)];
                    end
            end

            %                 h_la = figure; scatter(dataToFitToCircle(:, 1), dataToFitToCircle(:, 2)); hold on
            %                 lala = fit_ellipse( dataToFitToCircle(:, 1),dataToFitToCircle(:, 2), h_la );

            lala = fit_ellipse( dataToFitToCircle(:, 1),dataToFitToCircle(:, 2) );
            currPlateMean = [lala.X0 lala.Y0];
    end
end

function currPlateNorm = removeVarianceFct(obj, currInputUse)

    switch obj.removeVariance
        case 'none'
            currPlateNorm = ones(1, size(currInputUse, 2));

        case 'max'
            currPlateNorm = max(abs(currInputUse));

        case 'maxmax'
            currPlateNorm = max(max(abs(currInputUse))) * ones(1, size(currInputUse, 2));

        case 'circle'
            % offset is calculated from taking the circle fit
            switch obj.dataFormatName
                case 'PC'
                    dataToFitToCircle = currInputUse;
                    
                case 'PP'
                    switch obj.dataFormatDofStr
                        case {'all', 'alldof', 'sigdof'}
                            dataToFitToCircle = currInputUse;

                        case {'normalldof', 'normsigdof'}
                            % offset is calculated from taking the circle fit
                            dataToUse = currInputUse;
                            dataToFitToCircle = [dataToUse; dataToUse(:, 1) -dataToUse(:, 2)];
                    end
            end

            %                 h_la = figure; scatter(dataToFitToCircle(:, 1), dataToFitToCircle(:, 2)); hold on
            %                 lala = fit_ellipse( dataToFitToCircle(:, 1),dataToFitToCircle(:, 2), h_la );

            lala = fit_ellipse( dataToFitToCircle(:, 1),dataToFitToCircle(:, 2) );
            currPlateNorm = [lala.a lala.b];
    end
end

function [inputUse, currIndUse, dof] = sortDataPlates(obj, input, sigdofcellarray)
    switch obj.dataFormatName
        case 'PC'
            % now perform the dim reduction, since the DR may change from
            % instance to instance
            inputDR = obj.dimReduct.apply(input);
            currIndUse{1} = [];

            switch obj.dataFormatDofStr
                case 'all'
                    inputUse{1} = inputDR;
                    dof{1} = [0];

                case 'seq'
                    % doesn't apply for PCA
                    for i = 1:floor(size(inputDR, 2)/2)
                        inputUse{i} = inputDR(:, (i-1)*2+1:(i-1)*2+2);
                        dof{i} = [(i-1)*2+1:(i-1)*2+2];
                        currIndUse{i} = [];
                    end

                case 'fit'
                    %         if obj.dataFormatDof > 0
                    %             % restrict the DOFs
                    %             inputDR = inputDR(:, 1:obj.dataFormatDof);
                    %         end
            end

        case 'PP'
            allDofs = 1:size(input, 2)/2;

            % PP always only have 2 vectors (q and dq)
            inputQ  = input(:, allDofs);
            inputDQ = input(:, allDofs+max(allDofs));
            currIndUse{1} = [];

            switch obj.dataFormatDofStr
                case 'all' %all into one plate
                    inputUse{1} = input;
                    dof{1} = [allDofs];

                case 'alldof' % all into individual plates
                    %                 inputQ  = inputQ(:, :);
                    %                 inputDQ = inputDQ(:, :);

                    for i = 1:size(inputQ, 2)
                        inputUse{i} = [inputQ(:, i) inputDQ(:, i)];
                        dof{i} = [i];
                        currIndUse{i} = [];
                    end

                case 'sigdof'
                    %                 inputQ  = inputQ(:, obj.dataFormatDofNum);
                    %                 inputDQ = inputDQ(:, obj.dataFormatDofNum);

                    %invert the sigdofcellarray so each cell holds the ind for the
                    %data corresponding to that sigdof
                    if ~isempty(sigdofcellarray)
                        sigdofcellarrayInv = cell(length(allDofs), 1);
                        for i = 1:length(sigdofcellarray)
                            for j = 1:length(sigdofcellarray{i})
                                currSigDof = sigdofcellarray{i}(j);
                                sigdofcellarrayInv{currSigDof} = [sigdofcellarrayInv{currSigDof} i];
                            end
                        end
                    end

                    for i = 1:length(obj.dataFormatDofNum)
                        currSigDof = obj.dataFormatDofNum(i);
                        if ~isempty(sigdofcellarray)
                            currIndUse{i} = sigdofcellarrayInv{currSigDof};
                            inputUse{i} = [inputQ(currIndUse{i}, currSigDof) inputDQ(currIndUse{i}, currSigDof)];
                        else
                            % didn't pass in a sigdofcellarray so just use it
                            % all
                            currIndUse{i} = [];
                            inputUse{i} = [inputQ(:, currSigDof) inputDQ(:, currSigDof)];
                        end

                        dof{i} = [obj.dataFormatDofNum(i)];
                    end

                case 'normalldof' % of all the sigdofs
                    inputQ  = normVector(inputQ);
                    inputDQ = normVector(inputDQ);

                    inputUse{1} = [inputQ inputDQ];
                    dof{1} = [0];

                case 'normsigdof'
                    inputQ  = normVector(inputQ(:, obj.dataFormatDofNum));
                    inputDQ = normVector(inputDQ(:, obj.dataFormatDofNum));

                    inputUse{1} = [inputQ inputDQ];
                    dof{1} = [0];
            end
    end
    end