function [t, dist] = meanSVD_T(tCellOrig)
    
% %     dist = -1*ones(size(tCell, 2), size(tCell, 2));
    dist = [];
    tCell = checkTCell(tCellOrig);
    
    switch length(tCell)
        case 0
            % returns empty
            t = [];
            
        case 1
            % returns the passed in array
            t = tCell{1};
            
        otherwise
            % merge arrays
            [t, ind] = mergeTCell(tCell);
            
            R = t(1:3, 1:3);
            [S, V, D] = svd(R);
            V = eye(3);
            Rnew = S*V*D';
            
            T = t(1:3, 4);
            T = T / ind;
            
            t = eye(4);
            t(1:3, 1:3) = Rnew;
            t(1:3, 4) = T;
    end
end

function tCell = checkTCell(tCellOrig)
    % check the quality of the tCells
    tDist1 = [];
    tDist2 = [];
    tDist3 = [];
    tDistN = [];
    tCellTemp = {};
    for i = 1:length(tCellOrig)
        % remove empty tcells from consideration
        if isempty(tCellOrig{i}) || ~isempty(find(isnan(tCellOrig{i}), 1)) % if the cell is empty
            continue;
        end
        
         tCellTemp{end+1} = tCellOrig{i};
    end
    
    for i = 1:length(tCellTemp)
        % extra the relative frame distance individually
        if ~isempty(tCellTemp{i})
            tDist1 = [tDist1; tCellTemp{i}(1, 4)];
            tDist2 = [tDist2; tCellTemp{i}(2, 4)];
            tDist3 = [tDist3; tCellTemp{i}(3, 4)];
            tDistN = [tDistN; norm(tCellTemp{i}(1:3, 4))];
        end
    end

    [B1,TF1] = rmoutliers(tDist1);
    [B2,TF2] = rmoutliers(tDist2);
    [B3,TF3] = rmoutliers(tDist3);
    [BN,TFN] = rmoutliers(tDistN);
    
    meanTDist1 = mean(tDist1);
    meanTDist2 = mean(tDist2);
    meanTDist3 = mean(tDist3);
    meanTDistN = mean(tDistN);
    
    % in general, trust the individual ones
    tCell = {};
    baseThres = 0.03;
    for i = 1:length(tCellTemp)
        if TF1(i) == 0 && TF2(i) == 0 && TF3(i) == 0 %
            % if all 3 distances are not declared an outlier by rmoutliers'
            % MAD method
            tCell{end+1} = tCellTemp{i};
        elseif abs(tDist1(i) - meanTDist1) < baseThres && abs(tDist2(i) - meanTDist2) < baseThres && abs(tDist3(i) - meanTDist3) < baseThres
            % if the ranges tDist is small, the outlier rejection is very
            % narrow, so we can also trigger on a very small threshold
            tCell{end+1} = tCellTemp{i};
        elseif TFN(i) == 0 && abs(tDistN(i) - meanTDistN) < baseThres
            % if the norm is not declared an outlier, and within 3 cm of
            % the mean. the norm is generally a little less reliable
            % because the norm could be okay but the axes of the
            % individuals could be flipped
            tCell{end+1} = tCellTemp{i};
        end
    end
    
%     % however, if there is less than 2 individual distances, try the norm,
%     % but the norm might get stuck if there
%     if length(tCellTemp) > 2 && length(tCell) < 2
%         for i = 1:length(tCellTemp)
%             if TFN(i) == 0
%                 tCell{end+1} = tCellTemp{i};
%             end
%         end
%     end
end

function [t, ind] = mergeTCell(tCell)
    t = zeros(size(tCell{1}));
    ind = 0;
    
    for i = 1:length(tCell)
        if isempty(t) % if tcell{1} is empty
            t = zeros(size(tCell{i}));
        end

        ind = ind + 1;
        t = t + tCell{i};

        %         for j = 1:length(tCell)
        %             if isempty(tCell{j}) || ~isempty(find(isnan(tCell{j}), 1)) % if the cell is empty
        %                 continue;
        %             end
        %
        %             dist(i, j) = calcRotDiff(tCell{i}, tCell{j});
        %         end
    end
end