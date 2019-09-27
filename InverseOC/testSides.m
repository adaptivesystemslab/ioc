clc
startInd = 0;
for i = 1:10
%     % determine the current window to check
%     currFullWinIndsForward = startInd + (1:i);
%     
%     currFullWinIndsBackward = startInd - (1:i);
    
    
    addToRight = ceil((i-1)/2);
    addToLeft = ceil((i-2)/2);
    
    addInds = [-addToLeft:0 0:addToRight];
%     filtInd = addInds(addInds ~= 0);
    uniInds = unique(addInds);
    currFullWinIndsCenter = uniInds + startInd + 1
    
    
    
    
%     lala = [currFullWinIndsForward;
%         currFullWinIndsBackward;
%         currFullWinIndsCenter]    
end
