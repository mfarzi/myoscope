function addEchoTimes(obj, te)
    % addEchoTimes is a public method for class SCHEME
    
    assert(isempty(te)||isPositiveVector(te), 'MATLAB:scheme:addEchoTimes',...
        'Input argument must be a positive vector.');
    
    if isrow(te)
        te = te';
    end
    
    % remove nan values
    te(isnan(te)) = [];
    
    obj.teList = union(obj.teList, te, 'stable');
end
    