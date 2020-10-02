function select = pollDimReduction(name)

    select.settingName = name;
    parsing = regexp(name, '_', 'split');

    switch parsing{1}
        case 'None'
            select.name = 'None';

        case {'ksPCA', 'aksPCA'}
            select.name = 'ksPCA';
            sigmaVal = str2num(parsing{2})/10;

%             select.x_kernel = 'none';
%             select.y_kernel = 'delta';            
            select.sigma = sigmaVal;
            select.degree = str2num(parsing{3});            
            select.x_kernel = parsing{4};
            select.y_kernel = parsing{5};     
            select.LMatrixMultiplier = str2num(parsing{6});

            if select.degree == 0
                select.degree = -1;
            end
            
        case 'lasso'
            select.name = 'lasso';
            dimDeg = str2num(parsing{2});
            if dimDeg == 0
                dimDeg = -1;
            end
            select.degree = dimDeg;
            
        case 'Iso'
            select.name = 'Iso';
            
        case 'LMNN'
            select.name = 'LMNN';
            select.degree  = str2num(parsing{2});
            
        case {'LMNN2', 'aLMNN2'}
            select.name = 'LMNN2';
            select.degree  = str2num(parsing{2});
            
        case {'PCA', 'aPCA'}
            dimDeg = str2num(parsing{2});
            if dimDeg == 0
                dimDeg = -1;
            end
            
            select.name = 'PCA';
            select.degree = dimDeg;
            
        case 'fPCA'
            dimDeg = str2num(parsing{2});
            if dimDeg == 0
                dimDeg = -1;
            end
            
            select.name = 'fPCA';
            select.degree = dimDeg;
            
        case 'FDA'
             dimDeg = str2num(parsing{2});
            if dimDeg == 0
                dimDeg = -1;
            end
            
            select.name = 'FDA';
            select.degree = dimDeg;
    end
end