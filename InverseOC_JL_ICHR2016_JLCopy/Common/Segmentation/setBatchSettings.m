% batchSettings.dataSelectList = {'5_q-dq-ddq'};
batchSettings.dimReductionList = {'PCA_0', 'PCA_2', 'FDA_0', 'None'}; %, ...
%         'ksPCA_01_2', 'ksPCA_01_6', ...
%         'ksPCA_04_2', 'ksPCA_04_6', ...
%         'ksPCA_08_2', 'ksPCA_08_6', ...
%         'ksPCA_32_2', 'ksPCA_32_6', 'None'};
% batchSettings.dimReductionList = {'PCA_2', 'PCA_-1', 'FDA_-1', ...
%     'ksPCA_0.1_6', 'ksPCA_0.1_4', ...
%     'ksPCA_0.4_6', 'ksPCA_0.4_4', ...
%     'ksPCA_0.8_6', 'ksPCA_0.8_4'};
batchSettings.aggregatorList = {'None'...
        'Boosting_3', 'Bagging_3', ...
        'Boosting_5', 'Bagging_5'};
% batchSettings.classifierList = {'kNN_3', 'kNN_9', 'QDA', ...
%         'RBF_5', 'RBF_10', 'RBF_15', 'RBF_20', 'RBF_30', ...
%         'SVM_linear', 'SVM_polynomial', 'SVM_radial', 'SVM_sigmoid'...
%         'NN_10_10', 'NN_10_10_10', 'NN_20_10_5', 'NN_20_20_20'};
    
% batchSettings.classifierList = {'SVM_radial', 'SVM_polynomial', 'SVM_linear', 'SVM_sigmoid', ...
%         'RBF_5', 'RBF_10', 'RBF_15', 'RBF_20', 'RBF_30', ...
%         'kNN_3', 'kNN_9', 'QDA'};
    
batchSettings.classifierList = {'SVM_radial', 'SVM_polynomial', 'SVM_linear', 'SVM_sigmoid', ...
        'QDA'};