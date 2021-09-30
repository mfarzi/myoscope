function setNominalBvals(obj, b)
    % setNominalBvals is a public method for class SCHEME
    
    assert(isSemiPositiveVector(b), 'MATLAB:scheme:invalidInputArgument',...
        'Input argument must be a non-negative numeric vector.');
    
    if isrow(b)
        b = b';
    end
    
    obj.bvalDic = unique(b);
    [~,obj.bvalCode] = min(abs(obj.bval-obj.bvalDic'),[],2);
end
    