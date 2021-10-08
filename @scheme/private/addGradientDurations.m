function addGradientDurations(obj, delta)
    % addGradientDurations is a public method for class SCHEME
    
    assert(isPositiveVector(delta), 'MATLAB:scheme:addGradientDurations',...
        'Input argument must be a positive vector.');
    
    if isrow(delta)
        delta = delta';
    end
    
    obj.deltaList = union(obj.deltaList, delta, 'stable');
end
    