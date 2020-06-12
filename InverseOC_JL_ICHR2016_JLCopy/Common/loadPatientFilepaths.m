function [fileStack] = loadPatientFilepaths(sourcePath, specStruct)
    % loadPatientFilepaths compile a list of suitable patient filepaths,
    % based on the pathStruct specified
    
    blacklistPath = fullfile(sourcePath, '..', 'Blacklist.csv'); % if a file exists here, then use it

    if exist(blacklistPath, 'file') && isfield(specStruct, 'blackList') && specStruct.blackList
        blackListStruct = parseCSV_blacklist(blacklistPath);
    else
        blackListStruct = [];
    end
    
    fileStack = {};
    
    if ~exist('specStruct', 'var') || isempty(specStruct)
        % no specification. load all data available

    elseif isfield(specStruct, 'exercise')
        % legacy support, only 'exercise' is passed in
        if ~iscell(specStruct.exercise)
            % this function excepts the prefixes to be in cells
            tempCell{1} = specStruct.exercise;
        else
            tempCell = specStruct.exercise;
        end
        specStruct.exerciseAcceptPrefix = tempCell;
    end
    
    % make sure that missing specStruct fields are assigned blanks
    p = inputParser;
    p.KeepUnmatched = true;
    % replace the default values with the actual values if they're
    % being passed in explicitly
    addOptional(p, 'dataset', []);
    addOptional(p, 'patient', []);
    addOptional(p, 'session', []);
    addOptional(p, 'exerciseAcceptPrefix', []);
    addOptional(p, 'exerciseAcceptPrefixString', []);
    addOptional(p, 'exerciseAcceptPrefixInstance', []);
    addOptional(p, 'exerciseAcceptSuffix', []);
    addOptional(p, 'exerciseAcceptSuffixString', []);
    addOptional(p, 'exerciseAcceptSuffixInstance', []);
    addOptional(p, 'exerciseRejectPrefix', []);
    addOptional(p, 'exerciseRejectSuffix', []);
    addOptional(p, 'lastSetOnly', 0);
    
    parse(p,specStruct); % perform the file checking
    specStruct = p.Results;
    
    ind_fileStack = 0;
    
    dirPatient = dir(sourcePath);
    for ind_subj = 1:length(dirPatient)
        currSubj = dirPatient(ind_subj).name;

        if length(currSubj) < 7
            continue
        elseif ~strcmpi(currSubj(1:7), 'Subject')
            continue
        end

        currSubjNum = str2num(currSubj(8:end));

        currSubjPath = fullfile(sourcePath, currSubj);
        dirSession = dir(currSubjPath);

        for ind_sess = 1:length(dirSession)
            currSess = dirSession(ind_sess).name;

            if length(currSess) < 7
                continue
            elseif ~strcmpi(currSess(1:7), 'Session')
                continue
            end

            currSessNum = str2num(currSess(8:end));

            currSessPath = fullfile(currSubjPath, currSess);
            dirExercise = dir(currSessPath);
            
            for ind_exer = 1:length(dirExercise)
                currExer = dirExercise(ind_exer).name;

                if length(currExer) < 6
                    continue
                end
                
                if ~dirExercise(ind_exer).isdir
                    continue
                end
                
                % pull out the specific prefix/suffix for comparison
                % does the first instance of the constraints contain a
                % number? if it does, then we will need to check for the
                % instance number specifically. if not, then we will remove
                % it from the currExero
%                 if isempty(str2num(currExer(end)))
%                     % is the last digit a number?
%                     currExerType = currExer(1:end);
%                     currExerInst = 0;
%                 else
%                     currExerType = currExer(1:end-1);
%                     currExerInst = str2num(currExer(end));
%                 end
                
                currExerType = currExer(1:end);
                currExerInst = 0;
                
                % check to see if the instance is delimited by a underscore
                [currExerType, currExerInst] = exerciseStringDelimit(currExer);
                
                if specStruct.lastSetOnly 
                    % keep only the numerically last entry
                    otherExerType = dir(fullfile(currSessPath, [currExerType '*']));
                    
                    if ~strcmpi(currExer, otherExerType(end).name)
                        % it is not last one. keep going
                        continue
                    end
                end
                
                currExerStr = currExer(1:end);
                currExerPath = fullfile(currSessPath, currExer);

                passParam = checkParameters(specStruct, currSubjNum, currSessNum, currExer, currExerType, currExerInst);
                
                if ~passParam
                    continue 
                end
                
                % check for blacklist
                if ~isempty(blackListStruct)
                    blDataset = strcmpi(specStruct.dataset, blackListStruct.Dataset);
                    blSubject = strcmpi(num2str(currSubjNum), blackListStruct.Subject);
                    blSession = strcmpi(num2str(currSessNum), blackListStruct.Session);
                    blExercise = strcmpi(currExer, blackListStruct.Exercise);
                    
                    blTotal = blDataset + blSubject + blSession + blExercise;
                    
                    if max(blTotal) == 4 % has a match on all 4 fronts
                        fprintf('Blacklist entry detected: %s-Subj%u-Sess%u-%s. Entry discarded. \n', ... 
                            specStruct.dataset, currSubjNum, currSessNum, currExer);
                        continue
                    end
                end
                
                ind_fileStack = ind_fileStack + 1;
                fileStack{ind_fileStack}.subjectNumber = currSubjNum;
                fileStack{ind_fileStack}.sessionNumber = currSessNum;
                fileStack{ind_fileStack}.exerciseName = currExerStr;
                fileStack{ind_fileStack}.exerciseType = currExerType;
                fileStack{ind_fileStack}.filePath = currExerPath;
            end
        end
    end  
end

function passParam = checkParameters(specStruct, currSubjNum, currSessNum, currExerFullString, currExerType, currExerInst)
    passParam = 1;

    % apply px/sess/exercise specific filtering
    if ~isempty(specStruct.patient) && sum(currSubjNum == specStruct.patient) == 0
        % if there is an entry for 'specStruct.patient', but it is not
        % the current patient
        passParam = 0;
        return
    end
    
    if ~isempty(specStruct.session) 
        if ~iscell(specStruct.session) 
            if sum(currSessNum == specStruct.session) == 0
                % legacy usage. the session is not a cell, so we don't need to
                % search through it
                passParam = 0;
                return
            end
        else
            % the session var is a cell, so we need to find the proper
            % array. it should correspond to the patient interval
            currSubjNumInd = find(currSubjNum == specStruct.patient, 1);
            if ~isempty(currSubjNumInd) 
                if intersect(specStruct.session{currSubjNumInd}, currSessNum)
                % is okay
                else 
                    passParam = 0;
                    return
                end
            end
        end
    end
    
    % applies the the accept prefix condition
    if ~isempty(specStruct.exerciseAcceptPrefix)
        stringToCheck = checkLastDigit(specStruct.exerciseAcceptPrefix{1}, currExerType, currExerFullString);
        stringToCheck = stringToCheck(1:length(specStruct.exerciseAcceptPrefix{1}));
        
        if length(stringToCheck) >= length(specStruct.exerciseAcceptPrefix{1}) && ...
                ~sum(strcmpi(stringToCheck, specStruct.exerciseAcceptPrefix))
            passParam = 0;
        end
    end
    
    if ~isempty(specStruct.exerciseAcceptPrefixInstance)
        stringToCheck = currExerType(1:length(specStruct.exerciseAcceptPrefixString{1}));
        [stringCheckArray] = strcmp(stringToCheck, specStruct.exerciseAcceptPrefixString);
        
        if sum(stringCheckArray) == 0
            passParam = 0;
        else
            % may have passed
            stringCheckIndex = find(stringCheckArray);
            for ind = 1:length(stringCheckIndex)
                if specStruct.exerciseAcceptPrefixInstance(stringCheckIndex(ind)) == 0 || ...
                        specStruct.exerciseAcceptPrefixInstance(stringCheckIndex(ind)) == currExerInst
                    % passed
                    passParam = 1;
                    break
                else
                    passParam = 0;
                end
            end            
        end
    end

    % reject prefix
    if ~isempty(specStruct.exerciseRejectPrefix)
        stringToCheck = checkLastDigit(specStruct.exerciseRejectPrefix{1}, currExerNoInstance, currExerWithInstance);
        stringToCheck = stringToCheck(1:length(specStruct.exerciseRejectPrefix{1}));
        
        if sum(strcmp(stringToCheck, specStruct.exerciseRejectPrefix))
            passParam = 0;
        end
    end
    
    % accept suffix
    if ~isempty(specStruct.exerciseAcceptSuffix)
        suffixStringToCheck = checkLastDigit(specStruct.exerciseAcceptSuffix{1}, currExerNoInstance, currExerWithInstance);
        suffixStringToCheck = suffixStringToCheck(end-length(specStruct.exerciseAcceptSuffix{1}):end);
        
        if ~(sum(strcmp(suffixStringToCheck, specStruct.exerciseAcceptSuffix)))
            % applies the the accept suffix condition
            passParam = 0;
        end
    end
    
    if ~isempty(specStruct.exerciseAcceptSuffixInstance)
        stringToCheck = currExerType(end-length(specStruct.exerciseAcceptSuffixString{1})+1:end);
        stringCheckArray = strcmp(stringToCheck, specStruct.exerciseAcceptSuffixString);
        
        if sum(stringCheckArray) == 0
            passParam = 0;
        else
            % may have passed
            stringCheckIndex = find(stringCheckArray);
            for ind = 1:length(stringCheckIndex)
                if specStruct.exerciseAcceptSuffixInstance(stringCheckIndex(ind)) == 0 || ...
                        specStruct.exerciseAcceptSuffixInstance(stringCheckIndex(ind)) == currExerInst
                    % passed
                else
                    passParam = 0;
                end
            end
        end
    end
    
    % reject suffix
    if ~isempty(specStruct.exerciseRejectSuffix)
        suffixStringToCheck = checkLastDigit(specStruct.exerciseRejectSuffix{1}, currExerNoInstance, currExerWithInstance);
        suffixStringToCheck = suffixStringToCheck(end-length(specStruct.exerciseAcceptSuffix{1}):end);
        
        if sum(strcmp(suffixStringToCheck, specStruct.exerciseRejectSuffix))
            % applies the the accept suffix condition
            passParam = 0;
        end
    end
end

function stringToCheck = checkLastDigit(specStructString, currExerNoInstance, currExerWithInstance)

    if isempty(str2num(specStructString(end)))
        % not checking for instance number at the end
        checkInstance = 0;
        stringToCheck = currExerNoInstance;
    else
        % checking for instance number
        checkInstance = 1;

        % no instance
        stringToCheck = currExerWithInstance;
    end
end



% function checkOld
%     pass = 1;
% 
%     % apply px/sess/exercise specific filtering
%     if ~isempty(specStruct.patient) && sum(currSubjNum == specStruct.patient) == 0
%         % if there is an entry for 'specStruct.patient', but it is not
%         % the current patient
%         pass = 0;
%     elseif ~isempty(specStruct.session) && sum(currSessNum == specStruct.session) == 0
%         pass = 0;
%     elseif ~isempty(specStruct.exerciseAcceptPrefix) && length(currExer) >= length(specStruct.exerciseAcceptPrefix{1}) && ...
%             ~sum(strcmp(currExer(1:length(specStruct.exerciseAcceptPrefix{1})), specStruct.exerciseAcceptPrefix))
%         % applies the the accept prefix condition
%         pass = 0;
%     elseif ~isempty(specStruct.exerciseAcceptSuffix) && ...
%             ~(sum(strcmp(currExer(end-length(specStruct.exerciseAcceptSuffix{1}):end), specStruct.exerciseAcceptSuffix)) || ...
%             sum(strcmp(currExer(end-length(specStruct.exerciseAcceptSuffix{1}):end-1), specStruct.exerciseAcceptSuffix)))
%         % applies the the accept suffix condition
%         pass = 0;
%     elseif ~isempty(specStruct.exerciseRejectPrefix) && ...
%             sum(strcmp(currExer(1:length(specStruct.exerciseRejectPrefix{1})), specStruct.exerciseRejectPrefix))
%         % applies the the accept prefix condition
%         pass = 0;
%     elseif ~isempty(specStruct.exerciseRejectSuffix) && ...
%             (sum(strcmp(currExer(end-length(specStruct.exerciseRejectSuffix{1}):end), specStruct.exerciseRejectSuffix)) || ...
%             sum(strcmp(currExer(end-length(specStruct.exerciseRejectSuffix{1}):end-1), specStruct.exerciseRejectSuffix)))
%         % applies the the accept suffix condition
%         pass = 0;
%     end
% end
% 
