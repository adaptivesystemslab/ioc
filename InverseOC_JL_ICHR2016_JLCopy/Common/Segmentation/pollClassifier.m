function select = pollClassifier(name)

    select.settingName = name;
    parsing = regexp(name, '_', 'split');

    switch parsing{1}
        case 'PolyPCA'
            select.name = 'PolyPCA';
            
        case 'kNN'
            select.name = 'kNN';
            select.k = str2num(parsing{2});
            select.distanceMetric = 'Euclidean';

        case 'LDA'
            select.name = 'LDA';
            
        case 'QDA'
            select.name = 'QDA';
            
        case 'RBF'
            select.name = 'RBF';
            select.nodes = str2num(parsing{2});

        case 'Threshold'
            select.name = 'Threshold';
            select.percentage = str2num(parsing{2});
       
        case 'NN'
            select.name = 'NN';
            net_param = [];
            for i = 2:length(parsing)
                net_param = [net_param str2num(parsing{i})];
            end
            select.net_params = net_param;
            
        case 'SL'
            select.name = 'SL';
            
        case 'SVM'
            select.name = 'SVM';
            select.kernel = parsing{2};
            select.tuning = parsing{3};
            
        case 'SVMInc'
            select.name = 'SVMInc';
            
            select.kernel = parsing{2};
            select.tuning = parsing{3};
            select.resetRate = parsing{4};
            select.retrainRate = parsing{5};
            select.incLength = parsing{6};
            select.adaptMethod = parsing{7};
             
             adaptParam = {};
             for i = 8:length(parsing)
                 adaptParam{end+1} = parsing{i};
             end
             select.adaptParam = adaptParam;
            
        case 'SVMOneClass'
            select.name = 'SVMOneClass';
            select.kernel = parsing{2};  
            
        case 'SVMThreeClass'
            select.name = 'SVMThreeClass';
            select.kernel = parsing{2};  
            
        case 'SVMRetrain'
            select.name = 'SVMRetrain';
            select.kernel = parsing{2};
            select.tuning = parsing{3}; 
            
        case 'Compose'
            select.name = 'Compose';
            select.classifierName = parsing{2}; 
            
        case 'Cluster'
            select.name = 'Cluster';
            
            select.kernel = parsing{2};
            select.tuning = parsing{3};
            select.resetRate = parsing{4};
            select.retrainRate = parsing{5};
            select.incLength = parsing{6};
            select.adaptMethod = parsing{7};
            
            adaptParam = {};
            for i = 8:length(parsing)
                adaptParam{end+1} = parsing{i};
            end
            select.adaptParam = adaptParam;
    end
end