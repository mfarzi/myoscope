function updateGhatCode(obj)
    % updateGhatCode is a private method for class scheme
    % estimate ghatCode, i.e. closest direction from the codebook
    theta = acosd(abs(obj.ghat*obj.ghatDic'));
    [~, iMin] = min(theta,[],2);
    
    obj.ghatCode = iMin;
end