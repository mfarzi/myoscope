function setNominalDirections(obj, ghat)
    % setNominalDirections is a public method for class SCHEME
    
    assert(isDirectionMatrix(ghat), 'MATLAB:scheme:addNominalGhat',...
        'Input argument must be a matrix with 3 columns.');
    
    % normalise the norm to one
    ghatDic = scheme.normaliseDirections(unique(ghat,'row'));
    
    % if [0 0 0] exist in the ghatDic, it must be the first element.
    idx = find(ghatDic(:,1)==0 & ghatDic(:,2)==0 & ghatDic(:,3)==0);
    if idx>1
        ghatDic(idx,:) = [];
        ghatDic = [0 0 0;ghatDic];
    end
    obj.ghatDic = ghatDic;
    obj.updateGhatCode;
end
    